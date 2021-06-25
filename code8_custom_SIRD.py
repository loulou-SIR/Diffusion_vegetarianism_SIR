import networkx as nx
import matplotlib.pyplot as plt
import EoN
import xlrd
import random
from matplotlib import pyplot as plt
import numpy as np
import time
from collections import defaultdict
import heapq
import csv
from scipy.integrate import odeint


'''Custom code for SIRD stochastic model
1) Compartment D was added
2) The weights of the nodes were considered by modifying process_trans_SIR function
'''

##Classe myQueue

class myQueue(object):
    r'''
    This class is used to store and act on a priority queue of events for
    event-driven simulations.  It is based on heapq.
    Each queue is given a tmax (default is infinity) so that any event at later
    time is ignored.

    This is a priority queue of 4-tuples of the form
                   ``(t, counter, function, function_arguments)``
    The ``'counter'`` is present just to break ties, which generally only occur when
    multiple events are put in place for the initial condition, but could also
    occur in cases where events tend to happen at discrete times.
    note that the function is understood to have its first argument be t, and
    the tuple ``function_arguments`` does not include this first t.
    So function is called as
        ``function(t, *function_arguments)``
    Previously I used a class of events, but sorting using the __lt__ function
    I wrote was significantly slower than simply using tuples.
    '''
    def __init__(self, tmax=float("Inf")):
        self._Q_ = []
        self.tmax=tmax
        self.counter = 0 #tie-breaker for putting things in priority queue
    def add(self, time, function, args = ()):
        r'''time is the time of the event.  args are the arguments of the
        function not including the first argument which must be time'''
        if time<self.tmax:
            heapq.heappush(self._Q_, (time, self.counter, function, args))
            self.counter += 1
    def pop_and_run(self):
        r'''Pops the next event off the queue and performs the function'''
        t, counter, function, args = heapq.heappop(self._Q_)
        function(t, *args)
    def __len__(self):
        r'''this will allow us to use commands like ``while Q:`` '''
        return len(self._Q_)

##Function get rate functions

def _get_rate_functions_(G, tau, gamma1, gamma2, transmission_weight = None,
                        recovery_weight=None, death_weight=None): #add 2 custom parameters
    r'''
    Arguments :
        G : networkx Graph
            the graph disease spreads on
        tau : number
            disease parameter giving edge transmission rate (subject to edge scaling)
        gamma1: number (default None)
            disease parameter giving typical recovery rate,
        gamm2: death rate

        transmission_weight : string (default None)
            The attribute name under which transmission rates are saved.
            `G.adj[u][v][transmission_weight]` scales up or down the recovery rate.
            (note this is G.edge[u][v][..] in networkx 1.x and
            G.edges[u,v][..] in networkx 2.x.
            The backwards compatible version is G.adj[u][v]
            https://networkx.github.io/documentation/stable/release/migration_guide_from_1.x_to_2.0.html)
        recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their
            recovery rates
                `gamma_i = G.node[i][recovery_weight]*gamma`
    Returns :
        : trans_rate_fxn, rec_rate_fxn
            Two functions such that
            - `trans_rate_fxn(u,v)` is the transmission rate from u to v and
            - `rec_rate_fxn(u)` is the recovery rate of u.
'''
    if transmission_weight is None:
        trans_rate_fxn = lambda x, y: tau
    else:
        try:
            trans_rate_fxn = lambda x, y: tau*G.adj[x][y][transmission_weight] #x,y are nodes identifier
        except AttributeError: #apparently you have networkx v1.x not v2.x
            trans_rate_fxn = lambda x, y: tau*G.edge[x][y][transmission_weight]

    if recovery_weight is None:
        rec_rate_fxn = lambda x : gamma1
    else:
        rec_rate_fxn = lambda x : gamma1*G.nodes[x][recovery_weight]

    if death_weight is None:
        death_rate_fxn = lambda x : gamma2
    else:
        death_rate_fxn = lambda x : gamma2*G.nodes[x][death_weight]

    #custom line
    return trans_rate_fxn, rec_rate_fxn, death_rate_fxn

##Function truncated exponential, find trans and rec delays SIR, process trans SIR, process rec SIR, trans and rec time Markovian const trans

def _find_trans_and_rec_delays_SIR_(node, sus_neighbors, trans_time_fxn,
                                    rec_time_fxn, death_time_fxn, trans_time_args=(),
                                    rec_time_args=(), death_time_args=()):

    rec_delay = rec_time_fxn(node, *rec_time_args)
    death_delay = death_time_fxn(node, *death_time_args) #custom
    trans_delay={}
    for target in sus_neighbors:
        trans_delay[target] = trans_time_fxn(node, target, *trans_time_args)
    return trans_delay, rec_delay, death_delay #custom

def _truncated_exponential_(rate, T):
    r'''returns a number between 0 and T from an
    exponential distribution conditional on the outcome being between 0 and T'''
    t = random.expovariate(rate)
    L = int(t/T)
    return t - L*T

#This is the function we customized in the master's thesis
def _process_trans_SIR_custom(time, G, source, target, times, S, I, R, D, Q, status,
                            rec_time, death_time, pred_inf_time, transmissions,
                            trans_and_rec_time_fxn,
                            trans_and_rec_time_args = ()):
    r'''
    From figure A.4 of Kiss, Miller, & Simon.  Please cite the book if
    using this algorithm.
    :Arguments:
    time : number
        time of transmission
**G**  networkx Graph
    node : node
        node receiving transmission.
    times : list
        list of times at which events have happened
    S, I, R : lists
        lists of numbers of nodes of each status at each time
    Q : myQueue
        the queue of events
    status : dict
        dictionary giving status of each node
    rec_time : dict
        dictionary giving recovery time of each node
    death_time : dict
        dict giving death time of each node
    pred_inf_time : dict
        dictionary giving predicted infeciton time of nodes
    trans_and_rec_time_fxn : function
        trans_and_rec_time_fxn(node, susceptible_neighbors, *trans_and_rec_time_args)
        returns tuple consisting of
           dict of delays until transmission from node to neighbors and
           float having delay until recovery of node
        An example of how to use this appears in the code fast_SIR where
        depending on whether inputs are weighted, it constructs different
        versions of this function and then calls fast_nonMarkov_SIR.
    trans_and_rec_time_args : tuple (default empty)
        see trans_and_rec_time_fxn
    :Returns:

    nothing returned
    :MODIFIES:

    status : updates status of newly infected node
    rec_time : adds recovery time for node
    times : appends time of event
    S : appends new S (reduced by 1 from last)
    I : appends new I (increased by 1)
    R : appends new R (same as last)
    D : appends new D (same as last)
    Q : adds recovery and transmission events for newly infected node.
    pred_inf_time : updated for nodes that will receive transmission
    '''

    if status[target] == 'S': #nothing happens if already infected.
        status[target] = 'I'
        times.append(time)
        transmissions.append((time, source, target))
        S.append(S[-1]-1) #one less susceptible
        I.append(I[-1]+1) #one more infected
        R.append(R[-1])   #no change to recovered
        D.append(D[-1])  #no change to death #custom
        #Custom part of the code
        G.nodes[target]['vege']=1


        suscep_neighbors = [v for v in G.neighbors(target) if status[v]=='S']
        #Custom part of the code
        for v in suscep_neighbors:
            G[target][v]['weight']=weight_norm(G.nodes[v]['subcategory'])

        #custom
        trans_delay, rec_delay, death_delay= trans_and_rec_time_fxn(target, suscep_neighbors,
                                                *trans_and_rec_time_args)


        rec_time[target] = time + rec_delay
        death_time[target] = time + death_delay #custom
        '''
        if rec_time[target]<=Q.tmax:
            Q.add(rec_time[target], _process_rec_SIR_,
                            args = (target, times, S, I, R, status))
        '''
        #custom
        if rec_time[target]<death_time[target]:
            Q.add(rec_time[target], _process_rec_SIR_,
                            args = (target, times, S, I, R, D, status))
        else:
            Q.add(rec_time[target], _process_death_SIR_,
                            args = (target, times, S, I, R, D, status))

        #custom
        for v in trans_delay:
            inf_time = time + trans_delay[v]
            if inf_time<= rec_time[target] and inf_time<= death_time[target] and inf_time < pred_inf_time[v] and inf_time<=Q.tmax:
                Q.add(inf_time, _process_trans_SIR_custom,
                              args = (G, target, v, times, S, I, R, D, Q,
                                        status, rec_time, death_time, pred_inf_time,
                                        transmissions, trans_and_rec_time_fxn,
                                        trans_and_rec_time_args
                                     )
                             )
                pred_inf_time[v] = inf_time

#Custom
def _process_rec_SIR_(time, node, times, S, I, R, D, status):
    r'''From figure A.3 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    :Arguments:
        event : event
            has details on node and time
        times : list
            list of times at which events have happened
        S, I, R, D : lists
            lists of numbers of nodes of each status at each time
        status : dict
            dictionary giving status of each node
    :Returns:
        :
        Nothing
    MODIFIES
    ----------
    status : updates status of newly recovered node
    times : appends time of event
    S : appends new S (same as last)
    I : appends new I (decreased by 1)
    R : appends new R (increased by 1)
    D : appends new D (same as last)
    '''
    times.append(time)
    S.append(S[-1])   #no change to number susceptible
    I.append(I[-1]-1) #one less infected
    R.append(R[-1]+1) #one more recovered
    D.append(D[-1])   #no change to number of dead
    status[node] = 'R'
    #Custom part of the code
    #G.nodes[target]['vege']=0 #the node has recovered #not possible

#Custom
def _process_death_SIR_(time, node, times, S, I, R, D, status):
    r'''Custom function
    :Arguments:
        event : event
            has details on node and time
        times : list
            list of times at which events have happened
        S, I, R, D : lists
            lists of numbers of nodes of each status at each time
        status : dict
            dictionary giving status of each node
    :Returns:
        :
        Nothing
    MODIFIES
    ----------
    status : updates status of newly recovered node
    times : appends time of event
    S : appends new S (same as last)
    I : appends new I (decreased by 1)
    R : appends new R (same at least)
    D : appends new D (increased by 1)

    '''
    times.append(time)
    S.append(S[-1])   #no change to number susceptible
    I.append(I[-1]-1) #one less infected
    R.append(R[-1])   #no change to number recovered
    D.append(D[-1]+1) #increased by one
    status[node] = 'R'
    #Custom part of the code
    #G.nodes[target]['vege']=1 #the node will remain vegetarian now #not possible



def _trans_and_rec_time_Markovian_const_trans_(node, sus_neighbors, tau, rec_rate_fxn):
    r'''I introduced this with a goal of making the code run faster.  It looks
    like the fancy way of selecting the infectees and then choosing their
    infection times is slower than just cycling through, finding infection
    times and checking if that time is less than recovery time.  So I've
    commented out the more "sophisticated" approach.
    '''

    '''
    duration = random.expovariate(rec_rate_fxn(node)) #old line
    '''
    duration_rec = random.expovariate(rec_rate_fxn(node))
    duration_death= random.expovariate(death_rate_fxn(node))
    duration_min=min(duration_rec,duration_death)

    '''
    trans_prob = 1-np.exp(-tau*duration) #old line
    '''
    #trans_prob = 1-np.exp(-tau*duration_rec)*np.exp(-tau*duration_death) #custom new line
    trans_prob = 1-np.exp(-tau*duration_min)

    number_to_infect = np.random.binomial(len(sus_neighbors),trans_prob) #this line does a binomial game with the number of susceptible neighbors and the probability of infection
        #print(len(suscep_neighbors),number_to_infect,trans_prob, tau, duration)
    transmission_recipients = random.sample(sus_neighbors,number_to_infect)
    trans_delay = {}
    for v in transmission_recipients:
        #trans_delay[v] = _truncated_exponential_(tau, duration_rec+duration_death)#custom new line
        trans_delay[v] = _truncated_exponential_(tau, duration_min)
    return trans_delay, duration_min #idem
#     duration = random.expovariate(rec_rate_fxn(node))
#     trans_delay = {}
#
#
#     for v in sus_neighbors:
#         if tau == 0:
#             trans_delay[v] = float('Inf')
#         else:
#             trans_delay[v] = random.expovariate(tau)
# #        if delay<duration:
# #            trans_delay[v] = delay
#     return trans_delay, duration

##slow approach 1:
#    next_delay = random.expovariate(tau)
#    index, delay = int(next_delay//duration), next_delay%duration
#    while index<len(sus_neighbors):
#        trans_delay[sus_neighbors[index]] = delay
#        next_delay = random.expovariate(tau)
#        jump, delay = int(next_delay//duration), next_delay%duration
#        index += jump

##slow approach 2:
    #trans_prob = 1-np.exp(-tau*duration)
    #number_to_infect = np.random.binomial(len(sus_neighbors),trans_prob)
        #print(len(suscep_neighbors),number_to_infect,trans_prob, tau, duration)
    #transmission_recipients = random.sample(sus_neighbors,number_to_infect)
    #trans_delay = {}
    #for v in transmission_recipients:
    #    trans_delay[v] = _truncated_exponential_(tau, duration)
    return trans_delay, duration_min #idem


##Function fast non Markov SIR custom

def fast_nonMarkov_SIR_(G, trans_time_fxn=None,
                        rec_time_fxn=None, death_time_fxn=None,
                        trans_and_rec_time_fxn = None,
                        trans_time_args=(),
                        rec_time_args=(),
                        death_time_args=(),
                        trans_and_rec_time_args = (),
                        initial_infecteds = None,
                        initial_recovereds = None,
                        initial_deads = None,
                        rho=None, tmin = 0, tmax = float('Inf'),
                        return_full_data = False, sim_kwargs = None):
    r'''
    A modification of the algorithm in figure A.3 of Kiss, Miller, &
    Simon to allow for user-defined rules governing time of
    transmission.

    Please cite the book if using this algorithm.
    This is useful if the transmission rule is non-Markovian in time, or
    for more elaborate models.
    Allows the user to define functions (details below) to determine
    the rules of transmission times and recovery times.  There are two ways to do
    this.  The user can define a function that calculates the recovery time
    and another function that calculates the transmission time.  If recovery is after
    transmission, then transmission occurs.  We do this if the time to transmission
    is independent of the time to recovery.

    Alternately, the user may want to model a situation where time to transmission
    and time to recovery are not independent.  Then the user can define a single
    function (details below) that would determine both recovery and transmission times.

    :Arguments:
    **G** Networkx Graph

    **trans_time_fxn** a user-defined function
        returns the delay until transmission for an edge.  May depend
        on various arguments and need not be Markovian.

        Returns float
        Will be called using the form

        ``trans_delay = trans_time_fxn(source_node, target_node, *trans_time_args)``
            Here trans_time_args is a tuple of the additional
            arguments the functions needs.
        the source_node is the infected node
        the target_node is the node that may receive transmission
        rec_delay is the duration of source_node's infection, calculated
        by rec_time_fxn.
    **rec_time_fxn** a user-defined function
        returns the delay until recovery for a node.  May depend on various
        arguments and need not be Markovian.

        Returns float.
        Called using the form

        ``rec_delay = rec_time_fxn(node, *rec_time_args)``
            Here rec_time_args is a uple of additional arguments
            the function needs.

    **trans_and_rec_time_fxn** a user-defined function
        returns both a dict giving delay until transmissions for all edges
        from source to susceptible neighbors and a float giving delay until
        recovery of the source.

        Can only be used **INSTEAD OF** ``trans_time_fxn`` AND ``rec_time_fxn``.

        Gives an **ERROR** if these are also defined

        Called using the form
        ``trans_delay_dict, rec_delay = trans_and_rec_time_fxn(
                                           node, susceptible_neighbors,
                                           *trans_and_rec_time_args)``
        here trans_delay_dict is a dict whose keys are those neighbors
        who receive a transmission and rec_delay is a float.

    **trans_time_args** tuple
        see trans_time_fxn

    **rec_time_args** tuple
        see rec_time_fxn
    **trans_and_rec_time_args** tuple
        see trans_and_rec_time_fxn

    **initial_infecteds** node or iterable of nodes
        if a single node, then this node is initially infected

        if an iterable, then whole set is initially infected

        if None, then choose randomly based on rho.  If rho is also
        None, a random single node is chosen.

        If both initial_infecteds and rho are assigned, then there
        is an error.

    **initial_recovereds** iterable of nodes (default None)
        this whole collection is made recovered.

        Currently there is no test for consistency with initial_infecteds.

        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.

    **rho** number
        initial fraction infected. number is int(round(G.order()*rho))
    **tmin** number (default 0)
        starting time

    **tmax** number (default infinity)
        final time
    **return_full_data** boolean (default False)
        Tells whether a Simulation_Investigation object should be returned.
    **sim_kwargs** keyword arguments
        Any keyword arguments to be sent to the Simulation_Investigation object
        Only relevant if ``return_full_data=True``

    :Returns:

    **times, S, I, R** numpy arrays

    Or if ``return_full_data is True``

    **full_data**  Simulation_Investigation object
        from this we can extract the status history of all nodes
        We can also plot the network at given times
        and even create animations using class methods.

    :SAMPLE USE:

    ::


        import EoN
        import networkx as nx
        import matplotlib.pyplot as plt
        import random

        N=1000000
        G = nx.fast_gnp_random_graph(N, 5/(N-1.))


        #set up the code to handle constant transmission rate
        #with fixed recovery time.
        def trans_time_fxn(source, target, rate):
            return random.expovariate(rate)
        def rec_time_fxn(node,D):
            return D

        D = 5
        tau = 0.3
        initial_inf_count = 100
        t, S, I, R = EoN.fast_nonMarkov_SIR(G,
                                trans_time_fxn=trans_time_fxn,
                                rec_time_fxn=rec_time_fxn,
                                trans_time_args=(tau,),
                                rec_time_args=(D,),
                                initial_infecteds = range(initial_inf_count))

        # note the comma after ``tau`` and ``D``.  This is needed for python
        # to recognize these are tuples
        # initial condition has first 100 nodes in G infected.

    '''
    if rho and initial_infecteds:
        raise EoN.EoNError("cannot define both initial_infecteds and rho")
    if rho and initial_recovereds:
        raise EoN.EoNError("cannot define both initial_recovereds and rho")

    if (trans_time_fxn and not rec_time_fxn) or (rec_time_fxn and not trans_time_fxn):
        raise EoN.EoNError("must define both trans_time_fxn and rec_time_fxn or neither")
    elif trans_and_rec_time_fxn and trans_time_fxn:
        raise EoN.EoNError("cannot define trans_and_rec_time_fxn at the same time as trans_time_fxn and rec_time_fxn")
    elif not trans_and_rec_time_fxn and not trans_time_fxn:
        raise EoN.EoNError("if not defining trans_and_rec_time_fxn, must define trans_time_fxn and rec_time_fxn")

    if not trans_and_rec_time_fxn: #we define the joint function.
        trans_and_rec_time_fxn =  _find_trans_and_rec_delays_SIR_
        trans_and_rec_time_args = (trans_time_fxn, rec_time_fxn, death_time_fxn, trans_time_args, rec_time_args, death_time_args) #custom by adding third fxn and args

    #now we define the initial setup.
    status = defaultdict(lambda: 'S') #node status defaults to 'S'
    rec_time = defaultdict(lambda: tmin-1) #node recovery time defaults to -1
    death_time = defaultdict(lambda: tmin-1) #node death time defaults to -1 #custom
    if initial_recovereds is not None:
        for node in initial_recovereds:
            status[node] = 'R'
            rec_time[node] = tmin-1 #default value for these.  Ensures that the recovered nodes appear with a time
    if initial_deads is not None:
        for node in initial_deads:
            status[node] = 'D'
            death_time[node] = tmin-1 #default value for these.  Ensures that the recovered nodes appear with a time #custom
    pred_inf_time = defaultdict(lambda: float('Inf'))
        #infection time defaults to \infty  --- this could be set to tmax,
        #probably with a slight improvement to performance.

    Q = myQueue(tmax)

    if initial_infecteds is None:  #create initial infecteds list if not given
        if rho is None:
            initial_number = 1
        else:
            initial_number = int(round(G.order()*rho))
        initial_infecteds=random.sample(G.nodes(), initial_number)
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds]
    #else it is assumed to be a list of nodes.

    times, S, I, R, D= ([tmin], [G.order()], [0], [0], [0]) #custom by adding D and [0]
    transmissions = []

    for u in initial_infecteds:
        pred_inf_time[u] = tmin
        Q.add(tmin, _process_trans_SIR_custom, args=(G, None, u, times, S, I, R, D, Q,
                                                    status, rec_time, death_time,   #custom add death time
                                                    pred_inf_time, transmissions,
                                                    trans_and_rec_time_fxn,
                                                    trans_and_rec_time_args
                                                )
                        )

    #Note that when finally infected, pred_inf_time is correct
    #and rec_time is correct.
    #So if return_full_data is true, these are correct

    while Q:  #all the work is done in this while loop.
        Q.pop_and_run()

    #the initial infections were treated as ordinary infection events at
    #time 0.
    #So each initial infection added an entry at time 0 to lists.
    #We'd like to get rid these excess events.
    times = times[len(initial_infecteds):]
    S=S[len(initial_infecteds):]
    I=I[len(initial_infecteds):]
    R=R[len(initial_infecteds):]
    D=D[len(initial_infecteds):] #custom

    if not return_full_data: #we will often put ourselves in this situation
        return np.array(times), np.array(S), np.array(I), \
               np.array(R), np.array(D)
    else:
        #strip pred_inf_time and rec_time down to just the values for nodes
        #that became infected
        #could use iteritems for Python 2, by   try ... except AttributeError
        infection_times = {node:time for (node,time) in
                            pred_inf_time.items() if status[node]!='S'}
        recovery_times = {node:time for (node,time) in
                                rec_time.items() if status[node] =='R'}
        death_times = {node:time for (node,time) in
                                death_time.items() if status[node] =='D'} #custom


        node_history = _transform_to_node_history_(infection_times, recovery_times, death_times,
                                                    tmin, SIR = True)
        if sim_kwargs is None:
            sim_kwargs = {}
        return EoN.Simulation_Investigation(G, node_history, transmissions,
                                            possible_statuses = ['S', 'I', 'R', 'D'],
                                            **sim_kwargs)

##Function Fast SIR custom

def fast_SIR_custom2(G, tau, gamma1, gamma2, initial_infecteds = None, initial_recovereds = None, initial_deads = None,
                rho = None, tmin = 0, tmax=float('Inf'), transmission_weight = None,
                recovery_weight = None, death_weight = None, return_full_data = False, sim_kwargs = None):
    r'''
    fast SIR simulation for exponentially distributed infection and
    recovery times

    From figure A.3 of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.


    :Arguments:
    **G** networkx Graph
        The underlying network
    **tau** number
        transmission rate per edge
    **gamma** number
        recovery rate per node

    **initial_infecteds** node or iterable of nodes
        if a single node, then this node is initially infected

        if an iterable, then whole set is initially infected

        if None, then choose randomly based on rho.

        If rho is also None, a random single node is chosen.

        If both initial_infecteds and rho are assigned, then there
        is an error.

    **initial_recovereds** iterable of nodes (default None)
        this whole collection is made recovered.
        Currently there is no test for consistency with initial_infecteds.
        Understood that everyone who isn't infected or recovered initially
        is initially susceptible.
    **rho** number
        initial fraction infected. number is int(round(G.order()*rho))
    **tmin** number (default 0)
        starting time

    **tmax** number  (default Infinity)
        maximum time after which simulation will stop.
        the default of running to infinity is okay for SIR,
        but not for SIS.
    **transmission_weight**    string  (default None)
        the label for a weight given to the edges.
        transmission rate is
        G.adj[i][j][transmission_weight]*tau
    **recovery_weight**   string (default None))
        a label for a weight given to the nodes to scale their
        recovery rates
        gamma_i = G.nodes[i][recovery_weight]*gamma
    **return_full_data**   boolean (default False)
        Tells whether a Simulation_Investigation object should be returned.
    **sim_kwargs** keyword arguments
        Any keyword arguments to be sent to the Simulation_Investigation object
        Only relevant if ``return_full_data=True``

    :Returns:

    **times, S, I, R** numpy arrays

    Or if ``return_full_data is True``

    **full_data**  Simulation_Investigation object
            from this we can extract the status history of all nodes.
            We can also plot the network at given times
            and create animations using class methods.

    :SAMPLE USE:
    ::
        import networkx as nx
        import EoN
        import matplotlib.pyplot as plt

        G = nx.configuration_model([1,5,10]*100000)
        initial_size = 10000
        gamma = 1.
        tau = 0.3
        t, S, I, R = EoN.fast_SIR(G, tau, gamma,
                                    initial_infecteds = range(initial_size))

        plt.plot(t, I)
    '''
    #tested in test_SIR_dynamics
    if transmission_weight is not None or tau*gamma1*gamma2 == 0:
        trans_rate_fxn, rec_rate_fxn, death_rate_fxn = _get_rate_functions_(G, tau, gamma1, gamma2,
                                                    transmission_weight,
                                                    recovery_weight, death_weight)
        def trans_time_fxn(source, target, trans_rate_fxn):
            rate = trans_rate_fxn(source, target)
            if rate >0:
                return random.expovariate(rate)
            else:
                return float('Inf')
        def rec_time_fxn(node, rec_rate_fxn):
            rate = rec_rate_fxn(node)
            if rate >0:
                return random.expovariate(rate)
            else:
                return float('Inf')

        #custom
        def death_time_fxn(node, death_rate_fxn):
            rate = death_rate_fxn(node)
            if rate >0:
                return random.expovariate(rate)
            else:
                return float('Inf')

        trans_time_args = (trans_rate_fxn,)
        rec_time_args = (rec_rate_fxn,)
        death_time_args = (death_rate_fxn,)#custom
        return fast_nonMarkov_SIR_(G, trans_time_fxn = trans_time_fxn,
                        rec_time_fxn = rec_time_fxn,
                        death_time_fxn = death_time_fxn,
                        trans_time_args = trans_time_args,
                        rec_time_args = rec_time_args,
                        death_time_args = death_time_args,
                        initial_infecteds = initial_infecteds,
                        initial_recovereds = initial_recovereds,
                        initial_deads = initial_deads,
                        rho=rho, tmin = tmin, tmax = tmax,
                        return_full_data = return_full_data,
                        sim_kwargs=sim_kwargs)
    else:
        #the transmission rate is tau for all edges.  We can use this
        #to speed up the code.

        #get rec_rate_fxn (recovery rate may be variable)
        trans_rate_fxn, rec_rate_fxn, death_rate_fxn= _get_rate_functions_(G, tau, gamma1, gamma2,
                                                    transmission_weight,
                                                    recovery_weight, death_weight)

        return fast_nonMarkov_SIR_(G,
                        trans_and_rec_time_fxn=_trans_and_rec_time_Markovian_const_trans_,
                        trans_and_rec_time_args=(tau, rec_rate_fxn, death_rate_fxn),
                        initial_infecteds = initial_infecteds,
                        initial_recovereds = initial_recovereds,
                        initial_deads = initial_deads,
                        rho=rho, tmin = tmin, tmax = tmax,
                        return_full_data = return_full_data,
                        sim_kwargs=sim_kwargs)



##Function to produce many stochastic simulations and average them


def average_stochastic_custom(sol_number, beta, gamma1, gamma2):
    max_t=0
    max_len=0
    simulations=[]
    for i in range(sol_number):
        #at each simulation we generate a custom graph
        g=generate_custom(N1,K1,P1)
        #step that allows for initializing vegetarian nodes
        initial_vege=[v for v in g.nodes() if g.nodes[v]['vege']==1]

        #Use of our custom fast SIR function
        t1,S1,I1,R1, D1= fast_SIR_custom2(g, tau = beta, gamma1=gamma1, gamma2=gamma2, initial_infecteds=initial_vege, transmission_weight="weight", return_full_data=False)
        simulations.append([t1,S1,I1,R1,D1])

        if max(t1)>max_t:
            max_t=max(t1)
        if len(t1)>max_len:
            max_len=len(t1)

    #mean of the simulations
    t_mean,S_mean,I_mean,R_mean, D_mean,counter=np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000)
    for i in range(sol_number):
        sim=simulations[i]
        t,S,I,R,D=sim[0],sim[1],sim[2],sim[3],sim[4]

        new_t=np.linspace(min(t),max(t),5000)#creating a new time vector that will be fixed in number

        new_S=np.interp(new_t,t,S)#interpolating a new S vector from new_t, t, S
        new_I=np.interp(new_t,t,I)
        new_R=np.interp(new_t,t,R)
        new_D=np.interp(new_t,t,D)

        for j in range(5000):#calculating the sum
            t_mean[j]+=new_t[j]
            S_mean[j]+=new_S[j]
            I_mean[j]+=new_I[j]
            R_mean[j]+=new_R[j]
            D_mean[j]+=new_D[j]
            counter[j]+=1
    for k in range(5000): #calculating the mean value (replace by /5000)
        t_mean[k]/=counter[k]
        S_mean[k]/=counter[k]
        I_mean[k]/=counter[k]
        R_mean[k]/=counter[k]
        D_mean[k]/=counter[k]

    return t_mean,S_mean,I_mean,R_mean,D_mean


##Execute custom stochastic SIRD algorithm


#Network parameters
N1=10000
K1=4
P1=0.3
#SIR parameters
beta1=1
gamma1=0.1
gamma2=0.5
#Simulation number
sol_number=50

#Execute custom stochastic SIRD
t1,S1,I1,R1,D1= average_stochastic_custom(sol_number, beta1,gamma1,gamma2)


##Plot custom stochastic SIRD and deterministic SIRD on same figure

fig1 = plt.figure(facecolor='w')
fig1.suptitle('Custom SIRD algorithm N=10000 beta='+str(beta1)+' gamma1='+str(gamma1)+' gamma2='+str(gamma2))
ax1 = fig1.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#plot custom stochastic SIRD
ax1.plot(t1, S1/N1, 'b', linestyle='-', alpha=0.8, lw=2, label='Omnivorous')
ax1.plot(t1, I1/N1, 'r', linestyle='-', alpha=0.8, lw=2, label='Transitioning')
ax1.plot(t1, R1/N1, 'y', linestyle='-',alpha=0.8, lw=2, label='Ex-vegetarian')
ax1.plot(t1, D1/N1, 'g', linestyle='-',alpha=0.8, lw=2, label='Vegetarian')


ax1.set_xlabel('Time /days')
ax1.set_ylabel('Number (10000s)')
ax1.set_ylim(0,1.2)
ax1.yaxis.set_tick_params(length=0)
ax1.xaxis.set_tick_params(length=0)
ax1.grid(b=True, which='major', c='w', lw=2, ls='-')
legend1 = ax1.legend()
legend1.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax1.spines[spine].set_visible(False)

#show both figures
plt.show()

