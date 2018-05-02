Motion Profile Generators
====

Cubic
----

Built upon:

![q = a_0 + a_1t + a_2t^2 + a_3t^3](https://latex.codecogs.com/gif.latex?q%3Da_0&plus;a_1t&plus;a_2t%5E2&plus;a_3t%5E3)

With boundary conditions s.t. q = q0, v = v0 at t = t0 and q = qf, v = vf at t = tf. 

```cpp
motion::profile::cubic<float> MP;

// Motion profile from 0 to 5 in time 0 to 10
MP.set(0,10,0,5);

// Position at time 1
MP.q_at(1);

// Velocity at time 10
MP.v_at(10);

// Motion profile from last q (5) to 0 in time 10 to 15 
MP.set(15,0);
```
