# LatticeHopping
### **MC simulation of energy migration over β-NaYbF<sub>4</sub> lattice _via_ Dexter's exchange mechanism**

  The hexagonal-phase(β) NaREF<sub>4</sub> upconversion nanocrystals (RE = rare earth elements) are an important family of luminescent nanomaterials due to its unique optical properties. This piece of code simulates the energy migration over the Ytterbium sublattice network of NaYbF<sub>4</sub> _via_ Dexter's exchange mechanism. The β-NaYbF<sub>4</sub> crystal adopts a hexagonal lattice configuration where Yb3+ ions sites forms connected hexagons if viewed along the Z-direction, i.e. from the top of the unit cells. As you can see from the illustation of the unit cell, the solid green balls are Ytterbium ions (**1a** sites), half green/half orange balls are Sodium/Ytterbium ions with 50% occupancy each (**1f** sites) and half white/half orange balls are Sodium ions with 50% vacancies. When multiple unit cells were plotted together, we could easily visualize those connected hexagonal Yb<sup>3+</sup> networks. 
  
<img width="492" alt="image" src="https://user-images.githubusercontent.com/109502810/187412800-7f2070fa-57be-40c1-8211-a2a53169f6d1.png"><img width="486" alt="image" src="https://user-images.githubusercontent.com/109502810/187446251-d783c533-34d2-4f11-836c-a102bddba14c.png">


  In the Yb3+ sublattice, each Yb<sup>3+</sup> at the **1a** site has 5 cloest Yb<sup>3+</sup> neighbors (2 on **1a** sites and 3 on **1f** sites) while each Yb<sup>3+</sup> at the **1f** site has 6 cloest Yb<sup>3+</sup> neighbors (all on **1a** sites). The MC simulation starts with a quantum of photon energy absorbed by a Yb<sup>3+</sup> ion and migrates to its cloest neighbors with distance dependent probabilities described in [Dexter's original work](https://aip.scitation.org/doi/10.1063/1.1699044). Specifically speaking, the energy migrate from **1f** site to neighboring **1a** sites with equal probability of 1/6; the energy migrate from **1a** site to neighboring **1a** site with 44.35% and **1f** site with 3.76% probability.
  
<img width="500" alt="image" src="https://user-images.githubusercontent.com/109502810/187419941-8af0c997-5393-4464-8311-b6b0a87b4e9d.png">

We borrowed the [cube coordinates](https://www.redblobgames.com/grids/hexagons/) _(x,y,z)_ with _x+y+z=0_ to label the hexagons seen from the top. Locally each hexagon is not on the same plane but instead with 3 vertices on even integer height and 3 others on odd integer height in units of half of the unit cell height. Lattice grid points with odd integer height were the **1f** sites where the occupancy should be 50%. Therefore with (_x_, _y_, h) and the occupancy constaints, all the Yb3+ sites are well defined in space and are ready for MC simulation.

 <img width="700" alt="image" src="https://user-images.githubusercontent.com/109502810/187447845-51247804-0d3f-4216-8b4f-e9253c26a07a.png">

 
   The 
