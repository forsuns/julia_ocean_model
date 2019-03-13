# julia_ocean_model
The project that develop of a high performance ocean model using Julia language.
This ocean model based on Kämpf(2010)'s open source algorithms.

 In order to develop a high performance ocean model, we used Julia, a Just-In-Time compile language, and to obtain the solution of the momentum equation, we made the code to solve the Poisson equation by the Successive Over-Relaxation method. And the simplified Yellow Sea was modeled to test tidal forcing, Coriolis forces, and the effect of vertical eddy diffusivity coefficients. The model has been tested with a total of 5 cases in the scenario. As a result of the test, the model represented well the amphidromic points of Yellow Sea and the tidal characteristics of North West Sea and South West Sea of Korea. Therefore, it is considered that the ocean model using Julia language has developed successfully.
