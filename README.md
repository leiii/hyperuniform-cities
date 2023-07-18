# Hyperuniform organization in human settlements

Quantifying the spatial structure of human settlements (e.g., villages, towns, and cities) is fundamental to understanding the complexity of urban systems. However, it has long been unclear whether there is any hidden regularity in the spatial pattern of settlements that have emerged under the conditions of great geographic, cultural, and socioeconomic diversity. Here, using the geographic data from diverse regions, we find that the seemingly irregular distribution of settlements can be characterized by a simple geometric insight: hyperuniformity, an intriguing point configuration with small density fluctuations at large distances. Hyperuniformity has been identified in many physical and biological systems, but has rarely been documented in socioeconomic systems. To reproduce the settlement patterns, we develop a growth model with two intuitive mechanisms: spatial matching and competition. Our model accurately simulates how settlements evolve towards hyperuniformity and predicts the heavy-tail population distribution across settlements, in agreement with empirical observations. These results provide insights into the self-organization of cities, enrich the central place theory in human geography, and reveal the universality of spatial organization shared by social, physical, and biological systems.

***

Data:
- German village points (digitalized from Christaller's orginal central place map)
- China town points
- London MSOA points before and after the equal-density projection

Code:
- `hyperuniform_funcs.py` core functions to calculate density fluctuations, structure factors, and pair-correlation functions.

Paper:
- https://arxiv.org/abs/2306.04149
