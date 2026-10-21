### 01_thermal

This analytical model uses Fourier's law of heat conduction to model the delta across the potting material and heat sink base. It then uses the resulting surface temperature and the fin efficiency from a modified Bessel function. The Rayleigh and Nusselt numbers are then calculated based on the geometry to get the approximate convective/conductive heat transfer (Elenbaas). The radiated heat transfer is then calculated with the Stefan-Boltzmann law with a weighting factor for self-absorption to finally compute the total heat dissipation.
