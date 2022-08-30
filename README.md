# LatticeHopping
### **MC simulation of energy migration over β-NaYbF4 lattice _via_ Dexter's exchange mechanism**

  The hexagonal-phase(β) NaREF4 upconversion nanocrystals (RE = rare earth elements) are an important family of luminescent nanomaterials due to its unique optical properties. This piece of code simulates the energy migration over the Ytterbium sublattice network of NaYbF4 _via_ Dexter's exchange mechanism. The β-NaYbF4 crystal adopts a hexagonal lattice configuration where Yb3+ ions sites forms connected hexagons if viewed along the Z-direction, i.e. from the top of the unit cells. As you can see from the image of the unit cell shown here, the solid green balls are Ytterbium ions (**1a** sites), half green/half orange balls are Sodium/Ytterbium ions with 50% occupancy each (**1f** sites) and half white/half orange balls are Sodium ions with 50% vacancies. 
  
<img width="492" alt="image" src="https://user-images.githubusercontent.com/109502810/187412800-7f2070fa-57be-40c1-8211-a2a53169f6d1.png">

  In the Yb3+ sublattice, each Yb3+ at the **1a** site has 5 cloest Yb3+ neighbors (2 on **1a** sites and 3 on **1f** sites) while each Yb3+ at the **1f** site has 6 cloest Yb3+ neighbors (all on **1a** sites). The MC simulation starts with a quantum of photon energy absorbed by a Yb3+ ion and migrates to its cloest neighbors with distance dependent probabilities described in [Dexter's original work](https://aip.scitation.org/doi/10.1063/1.1699044). Specifically speaking, the energy migrate from **1f** site to neighboring **1a** sites with equal probability of 1/6; the energy migrate from **1a** site to neighboring **1a** site with 44.35% and **1f** site with 3.76% probability.
  
<img width="500" alt="image" src="https://user-images.githubusercontent.com/109502810/187419941-8af0c997-5393-4464-8311-b6b0a87b4e9d.png">

We borrowed the [cube coordinates](https://www.redblobgames.com/grids/hexagons/) _(x,y,z)_ with _x+y+z=0_ to label the hexagons seen from the top. Locally each hexagon is not on the same plane but instead with 3 vertices on even integer height and 3 others on odd integer height in units of half of the unit cell height. Lattice grid points with odd integer height were the **1f** sites where the occupancy should be 50%. Therefore with (_x_, _y_, h) and the occupancy constaints, all the Yb3+ sites are well defined in space and are ready for MC simulation.

 <img width="716" alt="image" src="https://user-images.githubusercontent.com/109502810/187438894-c66c7ea8-787d-4510-8602-e499cb5e9432.png">
 
   The 
