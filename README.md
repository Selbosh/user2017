useR!2017 poster
================
By David Selby and David Firth, Department of Statistics, University of Warwick.

A poster for the useR!2017 conference in Brussels.

Title
-----

Ranking influential communities in networks

Abstract
--------

Which scientific fields export the most intellectual influence, through recent research, to other fields? Citation behaviour varies greatly across disciplines, making inter-field comparisons difficult. Recent approaches to measuring influence in networks, based on the PageRank algorithm, take the source of citations (or recommendations or links) into account.

By aggregating all publications in each Web of Science field into "super-journals", we model the exchange of citations between fields. We then fit a Bradley--Terry paired comparisons model---which can be shown to be equivalent to scaled PageRank---to measure the relative influence of academic communities. Uncertainty intervals are provided to quantify the robustness of the ranking. All analyses were performed in R.

Which field is top of the table? Visit the poster to find out!

**Keywords**: PageRank, Bradley--Terry model, networks, ranking, bibliometrics

Design inspiration
------------------

Alberto Cairo's [two-page spread on the *Galileo* telescope](http://www.domusweb.it/content/dam/domusweb/en/interviews/2014/01/31/the_art_of_informationvisualization/rbig/Cairo2.jpg) is a good basis for a layout. This features one huge diagram just left of centre, with an introductory paragraph on the far left, a self-contained series of fun facts along the bottom and a set of more detailed diagrams on the right. See also [Giants of the Ocean](http://www.domusweb.it/content/dam/domusweb/en/interviews/2014/01/31/the_art_of_informationvisualization/rbig/Cairo1.jpg).

See also the diagram in Figure 3 of *Fast unfolding of communities in large networks* by [Blondel et al (2008)](https://arxiv.org/abs/0803.0476), where we have a large node-link diagram with communities aggregated into nodes, then a 'zoom' on a selected community, represented as its constituent nodes and links.

Colin Purrington's [results area design for portrait posters](https://i0.wp.com/colinpurrington.com/wp-content/uploads/2011/09/poster-template-vertical-2-purrington.jpg).

Method
======

We will need the following packages.

``` r
library(igraph)
library(ggraph)
```

Data wrangling
--------------

First we load the 2013 Web of Science citation data into R.

``` r
WoS <- readRDS('../thesis/data/thomsonreuters.Rds')
```

We are only interested in citations from 2013 to publications from 2003--2012. Furthermore, some miscellaneous journals are unhelpfully classified as `ALL OTHERS` and we want to omit any reference to these as well. Citations of count (`weight`) equal to zero should be deleted to avoid connectivity problems later.

``` r
WoS <- data.frame(from = WoS$Citing, to = WoS$Cited, weight = WoS$AllYears - WoS$Earlier)
WoS <- WoS[- union(grep('ALL OTHERS', WoS$from), grep('ALL OTHERS', WoS$to)), ]
WoS <- subset(WoS, weight > 0)
```

It's important to put the citing journal in the first column, the cited journal second and call the citation counts `weight`, because of how the `igraph` function `graph_from_data_frame` works. We will now turn our data into an `igraph` object.

``` r
ig <- graph_from_data_frame(WoS)
```

We now have a weighted, directed graph.

Before we run our community detection algorithm, there is some other housekeeping to do. A number of journals either gave out or received zero citations during the study period. These will break the graph up into disconnected components, so we want to remove all these singletons.

Firstly, we calculate the strongly-connected components of the graph. To be strongly connected to the rest of the network, a journal must both cite and be cited by others.

``` r
strong <- components(ig, mode = 'strong')
core <- which(strong$csize > 1)
if (length(core) > 1) stop('There should be only one core component')
```

The whole graph contains 12,604 nodes (journals). Of these, 10,713 belong to our "core" strongly-connected component and 1,891 are singletons, either weakly connected or completely disconnected from the rest of the graph. Let's get rid of these singleton nodes.

``` r
ig <- induced_subgraph(ig, which(strong$membership == core))
```

We now have 10,713 journals in our network.

Super journals
--------------

We will put our directed, weighted graph through the Infomap algorithm, as implemented in the `igraph` package. Results can be nondeterministic, so we will fix the random seed for reproducible results. Community detection can take a long time, so you might want to `cache` this chunk!

``` r
set.seed(2017)
infomap <- cluster_infomap(ig)
```

The algorithm returns 112 communities; the largest contains 971 journals and the smallest contains 2. The mean community size is 96 journals.

The journals of community will be aggregated into a *super-journal* representing all incoming and outgoing citations for that community. Edge weights are summed and other edge and vertex attributes are ignored.

``` r
sj <- contract.vertices(ig, membership(infomap), 'ignore')
sj <- simplify(sj, remove.multiple = TRUE, remove.loops = TRUE,
               edge.attr.comb = list(weight = 'sum', 'ignore'))
```

For later reference, each super-journal will be assigned a unique ID.

``` r
V(sj)$name <- 1:vcount(sj)
```

Visualisation
-------------

A nice way to visualise these graphs is with the `ggraph` package. Each node will be represented by a point, proportional to its PageRank score, and each edge will be represented by an arc between these points.

There is a problem, however: `ggraph`, `ggplot2` and our PDF reader won't like it if we try to plot 9,145,683 arcs in a single graphic!

(Another approach is to have a single arc for each pair of nodes, with opacity or width proportional to the number of citations, but this doesn't look very good in my opinion.)

Calculating the positions of nodes in our graphic will involve multidimensional scaling of the correlation matrix. Let's extract the weighted adjacency matrix, paying attention to the fact that `igraph` considers citations to travel from rows to columns.

``` r
xtab <- Matrix::t(as_adjacency_matrix(sj, attr = 'weight'))
```

We will scale the citations counts so no pair of nodes has more than 1000 arcs drawn between them on the final graphic.

``` r
small_xtab <- as.matrix(xtab)
diag(small_xtab) <- 0 # ignore self-citations
scalefactor <- 1000 / max(small_xtab)
small_xtab <- ceiling(small_xtab * scalefactor)
```

Then we can turn this back into an `igraph` object for visualisation. A useful attribute is the PageRank, which we can use to make nodes with greater total influence appear larger on the plot.

``` r
viz_ig <- graph_from_adjacency_matrix(t(small_xtab))
V(viz_ig)$PageRank <- page.rank(sj)$vector
```

Now let's perform multidimensional scaling to generate some coordinates.

``` r
my_mds <- create_layout(viz_ig,
                        layout = 'igraph',
                        algorithm = 'mds',
                        dist = 1 - cor(as.matrix(xtab)))
```

Let's make some plots!

``` r
extrafont::loadfonts(device = 'pdf', quiet = TRUE)

ggraph(my_mds) +
  geom_edge_fan0(alpha = .01, colour = 'steelblue3') +
  geom_node_point(aes(size = PageRank), fill = 'steelblue3', pch = 21, colour = 'white') +
  coord_fixed() +
  scale_x_reverse() + # flip horizontally
  theme_graph() +
  theme(legend.position = 'none')
```

On its own, our visualisation does not imply much because we don't know which community is which. We need labels on the communities for that. I have gone through and manually assigned plausible labels to the communities generated by Infomap.

``` r
labels <- read.csv('data/cluster_names.csv', stringsAsFactors = FALSE)
V(viz_ig)$field <- labels$field[match(V(viz_ig)$name, labels$community)]
```

We can then have a labelled plot for reference, or even show a selection of "interesting" labels on the main graph, while omitting most of them to avoid clutter.

``` r
extrafont::loadfonts(device = 'win', quiet = TRUE)

ggraph(my_mds) +
  geom_edge_fan0(alpha = .01, colour = 'steelblue3') +
  geom_node_point(aes(size = PageRank), fill = 'steelblue3', pch = 21, colour = 'white') +
  geom_node_text(aes(label = field), size = 3,
                 repel = TRUE,
                 family = 'Arial Narrow',
                 fontface = 'bold',
                 colour = 'steelblue',
                 segment.alpha = .2) +
  coord_fixed() +
  scale_x_reverse() +
  theme_graph() +
  theme(legend.position = 'none')
```

![](img/labelled-1.png)

Within-field analysis
---------------------

Now, let's take a particular field out of the network and examine its inner structure. Varin et al. (2016) set a precedent to study statistics journals, so let's have a look at statistics.

``` r
statistics <- labels$community[labels$field == 'statistics']
stats_ig <- induced_subgraph(ig, which(membership(infomap) == statistics))
V(stats_ig)$PageRank <- page.rank(stats_ig)$vector
```

There are 77 journals in the statistics subgraph, in 1 strongly-connected component(s).

We can visualise the network as before.

``` r
stats_xtab <- Matrix::t(as_adjacency_matrix(stats_ig))
stats_layout <- create_layout(stats_ig,
                              layout = 'igraph',
                              algorithm = 'mds',
                              dist = 1 - cor(as.matrix(stats_xtab)))
ggraph(stats_layout) +
  geom_edge_fan0(alpha = .05, colour = 'steelblue3') +
  geom_node_point(aes(size = PageRank), fill = 'steelblue3', pch = 21, colour = 'white') +
  geom_node_text(aes(label = name), size = 3,
                 repel = TRUE,
                 family = 'Arial Narrow',
                 fontface = 'bold',
                 colour = 'steelblue',
                 segment.alpha = .2) +
  coord_fixed() +
  theme_graph() +
  theme(legend.position = 'none')
```

![](img/stats_viz-1.png)
