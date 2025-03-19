# Figure explanation

# EDA 1 : general plots

## Fig1 

Bar plot. In x-axis is the number index assigned to the RNA family gathered from Rfam. Out of the 4300+ families, only the 329 with >50 sequences in the seed alignment are showed here.

What can be seen is that the overall gap frequency in the MSA is highly heterogenous, with however many families having alignments with over 50% gaps.

## Fig2

Showcases the standard deviation of gaps' frequency inside a family's MSA. most are around 5%, however a few manage to reach more than 10%.

## Fig3

The second figure shows a bar plot with the same x-axis as before. The y-axis is a gap count. To make this plot, we first consider the MSA of the family. In addition, we consider the consensus sequence for this family. For each sequence in the MSA, we then count the number of positions where the sequences possesses a nucleotide where a gap is expected (Hamming distance), i.e., the positions where a sequence **creates** a gap by being inherently different from the rest of the family. 

This count is averaged by the number of sequences in the family, giving the average number of gaps created per sequence in this family.

In the graph, we can see that for most families, only a few gaps are created by the sequences, however, some families appear prone to create lots of them, indicating the presence of disruptive sequences. Those sequences have the same function but very different sequence/structure than the rest of the family.

## Fig4

This graph is the continuity of figure 2 and shows the standard deviation of the number of gaps created per sequence in each family. x-axis is the same, y-axis is the standard deviation count of gaps created between sequences inside a single family.

Some families (also those with high gap creation per sequence) appear with deviations as high as their average, showcasing extreme variations between sequences. This concurs with the hypothesis of sequence's disruptiveness.


> The goal of the next graphs is to understand how isolated those disruptive sequences are. Do they tend to cause disruption alone or in small clusters with highly different sequence profile?

## Fig5

This graph isolates the "prominent fraction", i.e. the fraction of each family that deviates by more than $2\sigma$, $\sigma$ being the standard deviation. x-axis is the same, y-axis is the fraction.

Most families appear to have a fraction of about 5% above $2\sigma$, but some go as high as 17.5%. In addition, the families that possess the highest deviation in gaps per sequence do not appear to be those with the biggest prominent fraction. 

## Fig6 

Gap per sequence of the sequences in the prominent fraction. What is interesting here is that not all the family with a GPS much higher that the other have a prominent fraction with high GPS. Particularly, some families with normal GPS, showcase a much higher than average GPS in their prominent fraction (notably the far right families on the graph).

This suggests that in some families, it is indeed whole clusters of mildly different sequences (causing big prominent fraction but not an important number of GPS). On the contrary, the hypothesis of only a very small number of highly different sequences disrupting the alignment is supported as well (with small prominent fraction while high GPS in the prominent fraction).

## Fig7

Graph with the number of sequences in the alignment. The causation between high number of sequences and appearance of disruptive sequences is evident: more different sequences, means a higher probability of getting this kind of sequence.


## Fig8

The length of the alignment in the MSA does not appear to correlate, surprisingly. Intuitively, the bigger the alignment size, the more gaps, but this does not seem to be the case.

## Fig 9/10

The size variation of sequences seems to be typically around 10% on average, however, the maximum variation (i.e. between the biggest and smallest sequences of a family tipically goes from 20 to 60%)

## Fig11

Average size of sequences in each family.


# EDA 2 : correlations

## Fig1/2

Figures 1 and 2 display two correlations heatmap for each of the features in EDA 1. The first is focused on Pearson correlation, while the second is Spearman's. 

Figure 1 allows to see which features are strongly linearly connected. In this case we can see:
- Sequences' average size has high correlation with the GPS score creation and its deviation, as well as the length of the MSA (which is expected).
- The maximum variation in sequences' size is strongly correlated with the gaps frequency, the standard deviation of gaps' frequency and the deviation of sequences length (all three exoected). It is however not significantly correlated with GPS score of any kind. **Same conclusions are reached for the sequences' length variation**.
- The prominent fraction GPS score is strongly correlated with the number of sequences in the MSA. This is not unexpected : more sequences in a family implies more chances of disruptive sequences.
- The size of the prominent fraction does not correlated with annything significantly
- Deviation of the GPS score strongly correlates with sequences average size, just like the GPS score itself.

Take out message:
- The prominent fraction of GPS does NOT correlate with the sequence average size, while the GPS score and its deviation does. Thus, it is not sequences that creates the most gaps that allow the shortening of a sequence, but the rest. Thus, this proves that in the optics of shortening sequences, we could just ignore disruptive sequences if their size is not significantly smaller than that of the average msa size.

Figure 2 heatmap does not show differences in the correlations relationships.

## Fig3

The global gap fraction in the most significant alignments of Rfam is of about 30%.

## Fig4/5

The maximum length variation, gap frequency and deviation of the gap frequency do not correlate specifically with the sequences average size. However there is a strong linear correlation between the sequences length variation and the deviation in gaps frequency.


# EDA 3: gaps distribution

## Fig1

Figure 1 is a heatmap of the pairwise L2 distance norm between the gaps distributions frequency in Rfam families with more than 50 sequences in their MSA.

Gaps distribution are determined by counting the number of positions with $0, 1, 2, ..., n$ gaps for each families. Then each family has their 20-quantiles determined (i.e. each 5% of the distribution are averaged togerther). This allows to compare distributions between families despite them having different sequences size. The structure of the data is now a $N \times l$ matrix with l=20 and N=329, i.e. the number of families. Following this the pairwise matrix distance is computed as shown.

What can imediately be seen in this matrix is that it is relatively homogenous: most distributions have a distance of less than 0.6 with all of the others. However, a small portion of the families showcase strongly different distribution of their gaps respect to all others. They appear in red, and represent less than 1% of the sequences represented. Also, their pairwise distance is extremely low.

## Fig2

To confirm this, we check their gap distributions:
- They are all 12 identical
- present an absence of gaps with low frequency, and showcase exclusively gaps with very strong conservation, them being almost all in the last 5% quantile.

## Fig3 

A similar plotting is made for the rest of the heatmap (the most homogenous one). The outliers depicted below are taken out.

To increase the representativity of the depiction of the gap's distribution tendency for non-outlier families in Rfam, I computed the average distance for all families respect to all the other, and then partitioned the data in four groups according to this average in ascending order.

The line plot showcase thus the values of the average distribution of the group with the standard deviation. First of all, it is evident that the distributions are highly similar respect to one another, and moreover, respect to the 12 outliers. This confirms the intuition from the heatmap. In addition, the distribution shape hints at a possible separation in the distribution of gaps.

This hypothesis has to be prudent as this is only hinted by the curve's shape. However, there seems to be a "rupture" at some point in the curves, from whereon the slope abruptly increases. Hence, we hypothesize that this "rupture" could result from two types of gaps:
- The "indels", meaning gaps with sporadic appearances, caused by a non-deleterious insertion or deletion in other sequences
- The "chasm", meaning concomitant gaps appearing in a majority of sequences of the MSA due to a small number of sequences with highly divergent composition and/or structure.

The first quartile does not show any rupture: there is a possibility that the families in this quartile do not possess any chasm. The second, third and fourth show this rupture, and it could be hypothesized that the portion before it are made of indels, while after it are mostly chasm. 


Note that the granularity of the analysis is not optimal and that only a direct analysis family by family can determine whether there are or not chasm.

## Fig4

As a replacement of the GPS score, is the possibility of a sequence by sequence score. To do so, for each sequence in the MSA of a family, and for each position in those sequences where a nucleotide is present, we associate the gap frequency in the MSA at this position. Then, for each sequence this gap frequency is summed along the sequence's axis.

The idea is that if a sequence is very disruptive then, it will have a much higher score than the rest of the sequences in the alignment. Hereafter we computed this score for the outlier aforementioned.