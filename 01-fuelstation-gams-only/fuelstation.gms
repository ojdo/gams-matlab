$title Electric fuel station model (fuelstation.gms)

Sets
    t          time               / 1*8760 /
    i          type of production / pv, windon, windoff /
    tfirst(t)  first timestep
    tlast(t)   last timestep;

    tfirst(t) = yes$(ord(t) eq 1);
    tlast(t)  = yes$(ord(t) eq card(t));

Parameters
    cs      cost of storage tank (kEUR per MWh) / 100 /
    c(i)    cost of plant (kEUR per MW) / pv 3000, windon 1500, windoff 2500 /
    d(t)    demand (MWh)
    cf(t,i) relative (normalized to 1) production of plants;

    d(t) = uniform(0,1);
    cf(t,'pv') = min(max(0,power(sin(ord(t)/24*3.14/2),4)+normal(.1,.1)),1);
    cf(t,'windon') = min(max(0,uniform(0,1)),1);
    cf(t,'windoff')= min(max(0,sqrt(sqrt(uniform(0,1)))),1);

Variables
    x(i)       size of production facilities (MW)
    s          size of accumulator (MWh)
    st(t)      evolution of accumulator SOC (MWh)
    tp(t)      total production of plants per timestep (MWh)
    z          total cost (k�)
    sell(t);

Positive Variables x, s, st, sell;

Equations
    cost       total cost equation
    pp(t)      calculates tp (total production) from cf and x
    dd(t)      assures that demand is always satisfied
    storage(t) new_storage = storage + input - demand
    ss(t)      simulates the capacity of the accumulator
    ss0(t)     initial storage content
    ssN(t)     final storage content;

    cost..         z          =e= sum(i, x(i)*c(i)) + cs*s;
    pp(t)..        tp(t)      =e= sum(i, cf(t,i)*x(i));
    dd(t)..        d(t)       =l= tp(t) + st(t);
    storage(t)..   st(t+1)    =e= st(t) + tp(t) - d(t) - sell(t);
    ss(t)..        st(t)      =l= s;
    ss0(tfirst)..  st(tfirst) =g= s/2;
    ssN(tlast)..   st(tlast)  =g= s/2;

Model  fuelstation / all / ;
Solve  fuelstation using lp minimizing z ;
Display  x.l, s.l;
