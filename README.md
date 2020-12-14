# FuzzyQ
## Fuzzy Quantification of Common and Rare Species in Ecological Communities 
#### Juan A. Balbuena[&#x1f4e7;](mailto:j.a.balbuena@uv.es), Clara Monlleó-Borrull, Cristina Llopis-Belenguer, Isabel Blasco-Costa, Volodimir L. Sarabeev, Serge Morand
<p> Most species in ecological communities are rare whereas only a few are common. Although this distributional paradox has intrigued ecologists for decades, the interpretation of species abundance distributions remains elusive. 
Fuzzy Quantification of Common and Rare Species in Ecological Communities (FuzzyQ) shifts the focus from the prevailing species-categorization approach to develop a quantitative framework that seeks to place each species along a rare-commonness gradient. Given a community surveyed over a number of sites, quadrats, or any other convenient sampling unit, FuzzyQ uses a fuzzy clustering algorithm that estimates a probability for each species to be common or rare based on abundance-occupancy information (i.e., abundance and number of sites occupied). Such as probability can be interpreted as a commonness index ranging from 0 to 1. FuzzyQ also provides community-level metrics about the coherence of the allocation of species into the common and rare clusters that are informative of the nature of the community under study. </p>
<p> FuzzyQ produces ecological indicators easy to measure and interpret that can give both clear, actionable insights into the nature of ecological communities and provides a powerful way to monitor environmental change on ecosystems. Comparison among communities is greatly facilitated by the fact that the method is relatively independent of the number of sites or sampling units considered. FuzzyQ works satisfactorily with a wide range of communities varying in species richness, dispersion and abundance currencies. </p>
<p> The development version of FuzzyQ can be installed in R with the <code>devtools</code> package:</p>
<code>       
        library(devtools)<br>
        install_github("Ligophorus/FuzzyQ")<br>
</code>
<p> See also https://ligophorus.github.io/FuzzyQ/ </p>
