# julia_ocean_model
The project that develop of a high performance ocean model using Julia language

In order to develop a high performance ocean model, we used Julia, a Just-In-Time compile language, and to obtain the solution of the momentum equation, we made the code to solve the Poisson equation by the Successive Over-Relaxation method. And then we made two models to test Julia calculation codes.
First, a simple channel form is modeled to test constant source/sink conditions.
Second, the simplified Yellow Sea was modeled to test tidal forcing, Coriolis forces, and the effect of vertical eddy diffusivity coefficients.
The model has been tested with a total of eight cases in the two scenarios.
As a result of the test, the depth-averaged current speed of the three cases in Scenario 1 converged perfectly to the theoretical value, and that showed well a vertical flow velocity gradient due to the bottom friction.
Also, the result of Scenario 2 represented well the amphidromic points of Yellow Sea and the tidal characteristics of mid-western and southwestern coast of Korea.
Therefore, it is considered that the ocean model using Julia language has developed successfully, this suggests that the ocean model has come to the stage of successful transition from a classical compile language to a Just-In-Time compile language.
