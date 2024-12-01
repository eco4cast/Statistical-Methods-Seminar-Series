# Modeling Biological Processes as Stopped Random Walks

## Description

Many biological processes depend on the accumulation of a specific factor. For example, a flower first blooms in the spring after cumulative temperatures reach a threshold. (This is the widely used growing degree day model). A cumulative sum is a random walk, and thus these processes are examples of a stopped random walk. 

In this presentation:

- **Lizzie** will review examples from phenology where plant behavior is triggered by accumulations, such as forcing and chilling. See slides EFIstats_part1.pdf.
- **Jonathan** will provide a high-level overview of the Central Limit Theorem for stopped random walks. He will demonstrate how this large-sample approximation can be employed to model experimental and observational data in R. Additionally, Jonathan will discuss scenarios where the approximation may fail and introduce a non-asymptotic approach that offers a better fit. See directory [EFIstats_part2](https://github.com/eco4cast/Statistical-Methods-Seminar-Series/tree/main/auerbach-randomwalks/EFIstats_part2) containing slides, code, and data.

## Speakers

### **Dr. Jonathan Auerbach**
<img src="https://jauerbach.github.io/assets/images/me2.png" alt="Jonathan Auerbach" width="200">

**Statistician at George Mason University**

**Bio:**
Dr. Jonathan Auerbach is a statistician at George Mason University whose research covers a wide range of topics at the intersection of statistics and public policy. He has broad methodological interests in the analysis of longitudinal data, particularly for data science and causal inference. His policy interests include urban analytics, open data, and the collection, evaluation, and communication of official statistics.

### **Dr. Elizabeth Wolkovich**
<img src="https://biodiversity.ubc.ca/sites/default/files/styles/square_400/public/profile-images/lizzie.jpg" alt="Elizabeth Wolkovich" width="200">

**Ecologist at the University of British Columbia**

**Bio:**
Dr. Elizabeth Wolkovich is an ecologist at the University of British Columbia specializing in the ecological impacts of climate change on plant communities. She focuses on the impacts of warming on phenology—the timing of major growth and reproductive events—and its impacts on community assembly. Her research integrates ecological theory with cutting-edge statistical tools—often using Bayesian methods, including Stan (via RStan)—to predict how global climate shifts affect biodiversity and agriculture, particularly perennial crops like winegrapes.
