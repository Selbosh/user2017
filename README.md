useR!2017 poster
================
By David Selby and David Firth, Department of Statistics, University of Warwick.

A poster for the useR!2017 conference in Brussels.

### Title

Ranking influential communities in networks

### Abstract

Which scientific fields export the most intellectual influence, through recent research, to other fields? Citation behaviour varies greatly across disciplines, making inter-field comparisons difficult. Recent approaches to measuring influence in networks, based on the PageRank algorithm, take the source of citations (or recommendations or links) into account.

By aggregating all publications in each Web of Science field into "super-journals", we model the exchange of citations between fields. We then fit a Bradley–Terry paired comparisons model—which can be shown to be equivalent to scaled PageRank—to measure the relative influence of academic communities. Uncertainty intervals are provided to quantify the robustness of the ranking. All analyses were performed in R.

Which field is top of the table? Visit the poster to find out!

**Keywords**: PageRank, Bradley–Terry model, networks, ranking, bibliometrics

Design inspiration
------------------

Alberto Cairo's [two-page spread on the *Galileo* telescope](http://www.domusweb.it/content/dam/domusweb/en/interviews/2014/01/31/the_art_of_informationvisualization/rbig/Cairo2.jpg) is a good basis for a layout. This features one huge diagram just left of centre, with an introductory paragraph on the far left, a self-contained series of fun facts along the bottom and a set of more detailed diagrams on the right. See also [Giants of the Ocean](http://www.domusweb.it/content/dam/domusweb/en/interviews/2014/01/31/the_art_of_informationvisualization/rbig/Cairo1.jpg).

See also the diagram in Figure 3 of *Fast unfolding of communities in large networks* by [Blondel et al (2008)](https://arxiv.org/abs/0803.0476), where we have a large node-link diagram with communities aggregated into nodes, then a 'zoom' on a selected community, represented as its constituent nodes and links.

Since useR!2017 specifies posters must be A0 in portrait orientation, we can't implement landscape layouts like those mentioned above. While using some of those ideas, we might consider Colin Purrington's [results area design for portrait posters](https://i0.wp.com/colinpurrington.com/wp-content/uploads/2011/09/poster-template-vertical-2-purrington.jpg).

Another useful resource is the [Better Posters blog](http://betterposters.blogspot.co.uk/).

A large dose of inspiration may also be drawn from Dorling Kindersley *Eyewitness* books.

Method
======

Unfortunately, some of the raw data is property of Thomson Reuters / Clarivate Analytics so I cannot republish it in this repository. Nonetheless, hopefully you can get a general idea of what I have done from the source code and results.

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

We are only interested in citations from 2013 to publications from 2003–2012. Furthermore, some miscellaneous journals are unhelpfully classified as `ALL OTHERS` and we want to omit any reference to these as well. Citations of count (`weight`) equal to zero should be deleted to avoid connectivity problems later.

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
ggraph(my_mds) +
  geom_edge_fan0(alpha = .01, colour = '#4F94CD') +
  geom_node_point(aes(size = PageRank), fill = '#4F94CD', pch = 21, colour = 'white') +
  coord_fixed() +
  scale_x_reverse() + # flip horizontally
  theme_graph() +
  theme(legend.position = 'none')
```

On its own, our visualisation does not imply much because we don't know which community is which. We need labels on the communities for that. I have gone through and manually assigned plausible labels to the communities generated by Infomap.

``` r
labels <- read.csv('data/cluster_names.csv', stringsAsFactors = FALSE)
V(viz_ig)$field <- V(sj)$field <- labels$field[match(V(viz_ig)$name, labels$community)]
```

We can then have a labelled plot for reference, or even show a selection of "interesting" labels on the main graph, while omitting most of them to avoid clutter.

``` r
# Hide plot labels for very small fields
my_mds$flabel <- ifelse(rank(my_mds$PageRank) > 32, as.character(my_mds$field), NA)

ggraph(my_mds) +
  geom_edge_fan0(alpha = .007, colour = 'tomato2') + #4F94CD #2D7B95
  geom_node_point(aes(size = PageRank), fill = 'tomato2', pch = 21, colour = 'white') +
  geom_node_label(aes(label = flabel),
                  size = 2,
                  repel = TRUE,
                  family = 'Gill Sans MT Condensed',
                  fontface = 'bold',
                  colour = 'tomato2',
                  # Label options
                  segment.alpha = .5,
                  segment.size = 0.2,
                  label.r = unit(0.1, 'lines'),
                  label.size = NA, # no label border
                  label.padding = unit(0.1, 'lines'),
                  fill = rgb(1, 1, 1, .5)) +
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

At the moment we only have Thomson Reuters' abbreviations for the journal names. Some are more obvious than others! We have the full journal titles, but some are a bit too long, so we will use a table of custom short-ish titles.

``` r
titles <- read.csv('data/stats_titles.csv', stringsAsFactors = FALSE)
V(stats_ig)$title <- titles$Short[match(V(stats_ig)$name, titles$JCR)]
```

We can visualise the network as before.

``` r
stats_xtab <- Matrix::t(as_adjacency_matrix(stats_ig))
stats_layout <- create_layout(stats_ig,
                              layout = 'igraph',
                              algorithm = 'mds',
                              dist = 1 - cor(as.matrix(stats_xtab)))

ggraph(stats_layout) +
  geom_edge_fan0(alpha = .05, colour = 'tomato2') +
  geom_node_point(aes(size = PageRank), fill = 'tomato2', pch = 21, colour = 'white') +
  geom_node_text(aes(label = title), size = 3,
                 repel = TRUE,
                 family = 'Gill Sans MT Condensed',
                 fontface = 'bold',
                 colour = 'tomato2',
                 segment.alpha = .2) +
  coord_fixed() +
  theme_graph() +
  theme(legend.position = 'none')
```

![](img/stats_viz_web-1.png)

And here is a ranking of these statistical journals, by Bradley–Terry score.

``` r
(stats_ranks <- dplyr::arrange(
  data.frame(journal = V(stats_ig)$title,
             Scroogefactor = scrooge::Scroogefactor(stats_xtab),
             PageRank = V(stats_ig)$PageRank,
             BradleyTerry = scrooge::BTscores(stats_xtab),
             rank = rank(-scrooge::BTscores(stats_xtab), ties.method = 'first')),
  desc(BradleyTerry)
))
```

| journal                       |  Scroogefactor|   PageRank|  BradleyTerry|  rank|
|:------------------------------|--------------:|----------:|-------------:|-----:|
| JRSS-A                        |      0.0282376|  0.0103465|     0.0283302|     1|
| American Stat                 |      0.0270228|  0.0057930|     0.0270249|     2|
| JRSS-B                        |      0.0256999|  0.0454141|     0.0259910|     3|
| Canada J Stats                |      0.0239206|  0.0083131|     0.0241060|     4|
| Annals                        |      0.0236113|  0.1100324|     0.0240983|     5|
| Biometrika                    |      0.0228638|  0.0498885|     0.0230736|     6|
| JRSS-C                        |      0.0210658|  0.0089601|     0.0213421|     7|
| Technometrics                 |      0.0209600|  0.0139015|     0.0208967|     8|
| Stat Sci                      |      0.0209032|  0.0153774|     0.0206936|     9|
| Biostatistics                 |      0.0209463|  0.0158200|     0.0205330|    10|
| Lifetime Data Analysis        |      0.0182751|  0.0058858|     0.0199515|    11|
| J Stat Software               |      0.0201867|  0.0162305|     0.0197136|    12|
| JASA                          |      0.0194846|  0.0835435|     0.0194837|    13|
| Stats & Comp                  |      0.0199332|  0.0109201|     0.0193991|    14|
| JCGS                          |      0.0198372|  0.0190766|     0.0193987|    15|
| Scan J Stats                  |      0.0188394|  0.0169175|     0.0182769|    16|
| Biometrics                    |      0.0184227|  0.0371235|     0.0181814|    17|
| J Machine Learning Res        |      0.0178295|  0.0302504|     0.0172690|    18|
| Machine Learning              |      0.0164286|  0.0057000|     0.0171842|    19|
| Stat Meth Med Res             |      0.0168654|  0.0072661|     0.0167750|    20|
| Bernoulli                     |      0.0163526|  0.0165876|     0.0162539|    21|
| Int Stats Rev                 |      0.0161866|  0.0050933|     0.0162027|    22|
| Stat Neerlandica              |      0.0157788|  0.0038321|     0.0161538|    23|
| J Time Series                 |      0.0164238|  0.0075043|     0.0160605|    24|
| J Stat Plan & Infer           |      0.0166547|  0.0368382|     0.0159328|    25|
| ANZ J Stats                   |      0.0156540|  0.0046732|     0.0158607|    26|
| Stat Model                    |      0.0157907|  0.0047046|     0.0157993|    27|
| Ann Inst Stat Maths           |      0.0157902|  0.0083345|     0.0155305|    28|
| Extremes                      |      0.0149530|  0.0057294|     0.0153994|    29|
| Bayesian Analysis             |      0.0154541|  0.0114376|     0.0153939|    30|
| Test                          |      0.0162304|  0.0075401|     0.0152143|    31|
| J Qual Tech                   |      0.0133428|  0.0078417|     0.0144551|    32|
| Prob Eng & Inf Sci            |      0.0129198|  0.0038546|     0.0141473|    33|
| Biometrical J                 |      0.0149472|  0.0107017|     0.0141281|    34|
| Stats & Prob Letters          |      0.0146409|  0.0227675|     0.0137136|    35|
| Stats in Medicine             |      0.0145437|  0.0454388|     0.0136752|    36|
| J Nonpar Stats                |      0.0137078|  0.0080410|     0.0135831|    37|
| Comp Stats & Data Analysis    |      0.0146516|  0.0420382|     0.0134674|    38|
| J Multivariate Analysis       |      0.0142435|  0.0318751|     0.0132759|    39|
| Envir & Eco Stats             |      0.0129528|  0.0046764|     0.0131239|    40|
| Environometrics               |      0.0120375|  0.0072662|     0.0129543|    41|
| Scand Actuarial J             |      0.0119172|  0.0049733|     0.0124851|    42|
| Stoch Models                  |      0.0112702|  0.0041752|     0.0123314|    43|
| Ann Appl Stats                |      0.0135044|  0.0192067|     0.0120295|    44|
| Methods & Comp in Appl Prob   |      0.0108873|  0.0045467|     0.0115274|    45|
| Metrika                       |      0.0094645|  0.0055462|     0.0102697|    46|
| Stat Method                   |      0.0088005|  0.0046021|     0.0092631|    47|
| ESAIM Prob & Stats            |      0.0092176|  0.0036076|     0.0092150|    48|
| Comm Stats: Sim & Comp        |      0.0095515|  0.0060782|     0.0091717|    49|
| J Stats Comp & Sim            |      0.0092627|  0.0082753|     0.0091212|    50|
| Comm Stats: Theory & Methods  |      0.0094081|  0.0118277|     0.0090699|    51|
| Seq Analysis                  |      0.0088882|  0.0036688|     0.0089746|    52|
| EJ Stats                      |      0.0097919|  0.0154410|     0.0089544|    53|
| J Agr Bio & Envir Stats       |      0.0089445|  0.0045417|     0.0089017|    54|
| REVSTAT                       |      0.0086771|  0.0028658|     0.0083682|    55|
| SORT                          |      0.0079263|  0.0021687|     0.0083124|    56|
| Statistics                    |      0.0079254|  0.0050826|     0.0082635|    57|
| Stats & Interface             |      0.0087419|  0.0040690|     0.0079890|    58|
| J Applied Stats               |      0.0082399|  0.0066011|     0.0078485|    59|
| Appl Stoch Models in Business |      0.0067881|  0.0035229|     0.0077174|    60|
| Stats Papers                  |      0.0065067|  0.0060219|     0.0069310|    61|
| Qual Eng                      |      0.0063973|  0.0043052|     0.0068593|    62|
| Pakistan J Stats              |      0.0058241|  0.0040568|     0.0064576|    63|
| Int J Biostatistics           |      0.0073174|  0.0051220|     0.0063488|    64|
| Brazil J Prob & Stats         |      0.0056185|  0.0025517|     0.0059057|    65|
| Pharm Stats                   |      0.0055106|  0.0043995|     0.0058420|    66|
| J Biopharm Stats              |      0.0059426|  0.0077599|     0.0057170|    67|
| R Journal                     |      0.0063932|  0.0026105|     0.0054341|    68|
| Adv Stat Analysis             |      0.0051192|  0.0027805|     0.0052610|    69|
| Stat Methods & Appl           |      0.0044261|  0.0028058|     0.0050218|    70|
| Qual Rel Eng Int              |      0.0044372|  0.0060060|     0.0046339|    71|
| Comp Stats                    |      0.0048448|  0.0035324|     0.0046084|    72|
| J Korean SS                   |      0.0045534|  0.0032268|     0.0045248|    73|
| Stats in Biopharm Res         |      0.0038185|  0.0033177|     0.0044349|    74|
| Rev Colomb Estad              |      0.0041089|  0.0023813|     0.0043742|    75|
| Adv Data Analysis & Class     |      0.0032452|  0.0025590|     0.0036667|    76|
| Qual Tech & Quant Mgmt        |      0.0021391|  0.0022255|     0.0021176|    77|

And now for something completely different. What if we actually have a single non-statistics "super-journal" in the statistics network?

``` r
mapping <- match(V(ig)$name, V(stats_ig)$name)
mapping[is.na(mapping)] <- max(mapping, na.rm = TRUE) + 1
stats_others <- contract.vertices(ig, mapping, list(name = 'first', 'ignore'))
V(stats_others)$name[vcount(stats_others)] <- "(All others)"
others_xtab <- Matrix::t(as_adjacency_matrix(stats_others))
diag(others_xtab) <- 0

dplyr::arrange(
  data.frame(journal = V(stats_others)$name,
             Scroogefactor = scrooge::Scroogefactor(others_xtab),
             PageRank = page.rank(stats_others)$vector,
             BradleyTerry = scrooge::BTscores(others_xtab)),
  desc(BradleyTerry)
)
```

| journal              |  Scroogefactor|   PageRank|  BradleyTerry|
|:---------------------|--------------:|----------:|-------------:|
| J STAT SOFTW         |      0.0597587|  0.0159115|     0.0662567|
| J R STAT SOC B       |      0.0533165|  0.0320152|     0.0584122|
| AM STAT              |      0.0498796|  0.0084807|     0.0547457|
| BIOMETRIKA           |      0.0343575|  0.0368557|     0.0361016|
| J COMPUT GRAPH STAT  |      0.0275142|  0.0169518|     0.0283282|
| ANN STAT             |      0.0266800|  0.0694021|     0.0278383|
| STAT COMPUT          |      0.0268367|  0.0123934|     0.0275116|
| BIOMETRICS           |      0.0250801|  0.0277151|     0.0266368|
| J AM STAT ASSOC      |      0.0237224|  0.0560967|     0.0251152|
| MACH LEARN           |      0.0224210|  0.0070761|     0.0246029|
| J MACH LEARN RES     |      0.0212466|  0.0186240|     0.0225928|
| BIOSTATISTICS        |      0.0214237|  0.0136015|     0.0218223|
| CAN J STAT           |      0.0225656|  0.0100984|     0.0206261|
| STAT METHODS MED RES |      0.0188755|  0.0085729|     0.0203618|
| TECHNOMETRICS        |      0.0194398|  0.0148760|     0.0189859|
| STAT SCI             |      0.0190649|  0.0139816|     0.0189028|
| BIOMETRICAL J        |      0.0178737|  0.0118255|     0.0173150|
| J R STAT SOC C-APPL  |      0.0179341|  0.0102346|     0.0172444|
| J R STAT SOC A STAT  |      0.0171675|  0.0091808|     0.0170781|
| STAT MED             |      0.0158302|  0.0336120|     0.0167550|
| COMPUT STAT DATA AN  |      0.0161689|  0.0384759|     0.0164907|
| J QUAL TECHNOL       |      0.0141561|  0.0111754|     0.0154762|
| SCAND J STAT         |      0.0165599|  0.0157602|     0.0147462|
| ANN I STAT MATH      |      0.0154795|  0.0098950|     0.0141590|
| INT STAT REV         |      0.0142590|  0.0073603|     0.0136896|
| J STAT PLAN INFER    |      0.0141285|  0.0328774|     0.0133039|
| EXTREMES             |      0.0133718|  0.0081363|     0.0128179|
| LIFETIME DATA ANAL   |      0.0127270|  0.0081576|     0.0125918|
| BAYESIAN ANAL        |      0.0137490|  0.0115487|     0.0124839|
| J TIME SER ANAL      |      0.0137731|  0.0093792|     0.0124744|
| STAT MODEL           |      0.0132021|  0.0074011|     0.0117373|
| TEST-SPAIN           |      0.0128800|  0.0098124|     0.0102774|
| J BIOPHARM STAT      |      0.0096748|  0.0098427|     0.0101378|
| J NONPARAMETR STAT   |      0.0123958|  0.0098934|     0.0098247|
| ENVIRONMETRICS       |      0.0096255|  0.0093404|     0.0097203|
| AUST NZ J STAT       |      0.0108198|  0.0074397|     0.0093956|
| PHARM STAT           |      0.0083334|  0.0073263|     0.0093387|
| STAT NEERL           |      0.0104667|  0.0067418|     0.0090792|
| ENVIRON ECOL STAT    |      0.0088240|  0.0071893|     0.0089873|
| PROBAB ENG INFORM SC |      0.0081377|  0.0066973|     0.0086932|
| J MULTIVARIATE ANAL  |      0.0097456|  0.0262258|     0.0086846|
| SORT-STAT OPER RES T |      0.0074250|  0.0054547|     0.0081761|
| REVSTAT-STAT J       |      0.0079691|  0.0061630|     0.0079729|
| METRIKA              |      0.0077236|  0.0085714|     0.0078223|
| BERNOULLI            |      0.0086403|  0.0148244|     0.0073900|
| J AGR BIOL ENVIR ST  |      0.0076118|  0.0071660|     0.0072628|
| SCAND ACTUAR J       |      0.0068412|  0.0073625|     0.0072348|
| STAT PROBABIL LETT   |      0.0082754|  0.0212788|     0.0071910|
| METHODOL COMPUT APPL |      0.0066574|  0.0071288|     0.0069508|
| STOCH MODELS         |      0.0061108|  0.0066631|     0.0065910|
| STATISTICS           |      0.0068530|  0.0081674|     0.0064932|
| STAT METHODOL        |      0.0068281|  0.0076906|     0.0063508|
| (All others)         |      0.0056965|  0.0115052|     0.0061721|
| QUAL ENG             |      0.0057950|  0.0076824|     0.0061600|
| COMMUN STAT-SIMUL C  |      0.0071682|  0.0090134|     0.0060598|
| J STAT COMPUT SIM    |      0.0063426|  0.0112566|     0.0058592|
| SEQUENTIAL ANAL      |      0.0067499|  0.0071786|     0.0057831|
| COMMUN STAT-THEOR M  |      0.0059140|  0.0140319|     0.0054953|
| R J                  |      0.0065317|  0.0059580|     0.0052631|
| QUAL RELIAB ENG INT  |      0.0049446|  0.0092773|     0.0052533|
| INT J BIOSTAT        |      0.0062151|  0.0072493|     0.0052001|
| REV COLOMB ESTAD     |      0.0039183|  0.0060033|     0.0045831|
| ANN APPL STAT        |      0.0057245|  0.0152338|     0.0045302|
| ESAIM-PROBAB STAT    |      0.0048405|  0.0063724|     0.0045107|
| STAT BIOPHARM RES    |      0.0040342|  0.0063369|     0.0043468|
| STAT PAP             |      0.0041720|  0.0095298|     0.0040877|
| BRAZ J PROBAB STAT   |      0.0039294|  0.0057803|     0.0039361|
| STAT INTERFACE       |      0.0054177|  0.0067122|     0.0038533|
| APPL STOCH MODEL BUS |      0.0035906|  0.0064436|     0.0038035|
| ASTA-ADV STAT ANAL   |      0.0043451|  0.0061006|     0.0037645|
| J APPL STAT          |      0.0038169|  0.0094288|     0.0035409|
| ADV DATA ANAL CLASSI |      0.0035820|  0.0058913|     0.0031111|
| ELECTRON J STAT      |      0.0049902|  0.0129883|     0.0030936|
| QUAL TECHNOL QUANT M |      0.0030537|  0.0055776|     0.0026304|
| STAT METHOD APPL-GER |      0.0025070|  0.0061124|     0.0025973|
| J KOREAN STAT SOC    |      0.0024670|  0.0063257|     0.0021873|
| COMPUTATION STAT     |      0.0025795|  0.0064853|     0.0021572|
| PAK J STAT           |      0.0012715|  0.0068723|     0.0012420|

We can try plotting it, too.

``` r
stats_others <- simplify(stats_others, remove.loops = TRUE, remove.multiple = FALSE) # remove self-citations
V(stats_others)$PageRank <- page.rank(stats_others)$vector

other_layout <- create_layout(stats_others,
                              layout = 'igraph',
                              algorithm = 'mds',
                              dist = 1 - cor(as.matrix(others_xtab)))

ggraph(other_layout) +
  geom_edge_fan0(alpha = .01, colour = '#4F94CD') +
  geom_node_point(aes(size = PageRank), fill = '#4F94CD', pch = 21, colour = 'white') +
  geom_node_text(aes(label = name), size = 3,
                 repel = TRUE,
                 family = 'Gill Sans MT Condensed',
                 fontface = 'bold',
                 colour = '#4F94CD',
                 segment.alpha = .2) +
  coord_fixed() +
  theme_graph() +
  theme(legend.position = 'none')
```

![](img/plot_stats_others-1.png)

Ranking
-------

We can penalise the Bradley–Terry model by adding a "player zero" who cites/is cited by every player/journal/field at a constant rate (say 1/2). This will help reduce the chance of outliers (such as journals for which we have very little citation data, or fields containing very few journals) from shooting to the top or the bottom of a Bradley–Terry scores league table.

``` r
zero_cite <- .15 * sum(xtab) / 2 / nrow(xtab) ## 2ka = .15n
penalised_xtab <- rbind(zero_cite, cbind(zero_cite, xtab)) # add zeroth player

library(scrooge)
field_ranks <- data.frame(
  field = V(sj)$field,
  PageRank = PageRank(penalised_xtab)[-1],
  BradleyTerry = BTscores(penalised_xtab)[-1],
  Scroogefactor = Scroogefactor(penalised_xtab)[-1]
)
field_ranks$rank <- rank(-field_ranks$BradleyTerry, ties.method = 'first')
dplyr::arrange(field_ranks, desc(BradleyTerry))
```

| field                             |   PageRank|  BradleyTerry|  Scroogefactor|  rank|
|:----------------------------------|----------:|-------------:|--------------:|-----:|
| economics                         |  0.0095262|     0.0198675|      0.0198085|     1|
| statistics                        |  0.0052172|     0.0182253|      0.0184844|     2|
| rheumatology                      |  0.0075946|     0.0154524|      0.0149113|     3|
| medicine                          |  0.0785679|     0.0147539|      0.0140610|     4|
| surgery                           |  0.0187227|     0.0147393|      0.0140006|     5|
| nephrology                        |  0.0091929|     0.0147153|      0.0141108|     6|
| psychometrics                     |  0.0030642|     0.0146652|      0.0153434|     7|
| oncology                          |  0.0298477|     0.0142217|      0.0134646|     8|
| neuroscience                      |  0.0483922|     0.0138453|      0.0131796|     9|
| urology                           |  0.0054990|     0.0134027|      0.0129622|    10|
| politics                          |  0.0050798|     0.0129090|      0.0134032|    11|
| orthopaedics                      |  0.0081823|     0.0127964|      0.0125031|    12|
| parasitology                      |  0.0200139|     0.0127915|      0.0122454|    13|
| biomedical sciences               |  0.1391158|     0.0126568|      0.0119706|    14|
| psychiatry                        |  0.0162007|     0.0125158|      0.0119816|    15|
| psychology                        |  0.0157250|     0.0123666|      0.0119232|    16|
| management                        |  0.0082283|     0.0123594|      0.0118723|    17|
| vascular surgery                  |  0.0032286|     0.0120399|      0.0120587|    18|
| obstetrics                        |  0.0077170|     0.0117973|      0.0113524|    19|
| opthalmology                      |  0.0049835|     0.0117134|      0.0113681|    20|
| sports medicine                   |  0.0062347|     0.0115786|      0.0112158|    21|
| marketing                         |  0.0044868|     0.0112366|      0.0109821|    22|
| sociology                         |  0.0065198|     0.0109829|      0.0108478|    23|
| computer graphics                 |  0.0028843|     0.0106927|      0.0107898|    24|
| plastic surgery                   |  0.0036232|     0.0104481|      0.0102228|    25|
| radiology                         |  0.0104202|     0.0101617|      0.0098458|    26|
| dentistry                         |  0.0048420|     0.0095513|      0.0092427|    27|
| operational research              |  0.0058023|     0.0094057|      0.0092595|    28|
| rhinology                         |  0.0048558|     0.0093847|      0.0091205|    29|
| philosophy                        |  0.0024042|     0.0093150|      0.0098626|    30|
| zoology                           |  0.0212266|     0.0093070|      0.0090279|    31|
| communication                     |  0.0032187|     0.0090557|      0.0091353|    32|
| dermatology                       |  0.0046755|     0.0089293|      0.0085348|    33|
| social anthropology               |  0.0027730|     0.0088286|      0.0090860|    34|
| mathematical finance              |  0.0023882|     0.0087307|      0.0091967|    35|
| music                             |  0.0023110|     0.0086923|      0.0090142|    36|
| hypnosis                          |  0.0022330|     0.0086700|      0.0091130|    37|
| history of science                |  0.0025406|     0.0086569|      0.0088251|    38|
| evaluation                        |  0.0022738|     0.0086527|      0.0090182|    39|
| toxicology                        |  0.0082432|     0.0086072|      0.0080802|    40|
| lighting                          |  0.0022472|     0.0085945|      0.0089657|    41|
| psychoanalysis                    |  0.0022574|     0.0085334|      0.0088749|    42|
| transfusion                       |  0.0027156|     0.0085100|      0.0084292|    43|
| social history                    |  0.0025049|     0.0084987|      0.0088589|    44|
| French sociology                  |  0.0022331|     0.0084802|      0.0088812|    45|
| rehabilitation                    |  0.0022275|     0.0084746|      0.0088672|    46|
| onomastics/history of mathematics |  0.0022206|     0.0084621|      0.0088660|    47|
| natural history                   |  0.0022199|     0.0084576|      0.0088621|    48|
| Russian sociology                 |  0.0022237|     0.0084546|      0.0088548|    49|
| history of education              |  0.0022232|     0.0084451|      0.0088458|    50|
| history of economics              |  0.0022352|     0.0084432|      0.0088541|    51|
| religion                          |  0.0022335|     0.0084329|      0.0087993|    52|
| complexity                        |  0.0031541|     0.0084078|      0.0086762|    53|
| translation                       |  0.0022264|     0.0083900|      0.0087861|    54|
| agriculture                       |  0.0139732|     0.0083396|      0.0081027|    55|
| German sociology                  |  0.0022213|     0.0083082|      0.0086561|    56|
| linguistics                       |  0.0027522|     0.0082932|      0.0083216|    57|
| civil engineering                 |  0.0022271|     0.0082824|      0.0086966|    58|
| petrochemical engineering         |  0.0023708|     0.0082457|      0.0085083|    59|
| Brazilian/French philosophy       |  0.0022209|     0.0082120|      0.0084957|    60|
| leisure studies                   |  0.0023499|     0.0082049|      0.0084187|    61|
| Latin American studies            |  0.0022265|     0.0081671|      0.0085518|    62|
| Romanian/Italian ethics           |  0.0022224|     0.0081632|      0.0084943|    63|
| leather                           |  0.0022236|     0.0081602|      0.0085846|    64|
| anthropology                      |  0.0032005|     0.0081559|      0.0080893|    65|
| policy                            |  0.0027347|     0.0080814|      0.0081088|    66|
| international law                 |  0.0024032|     0.0080547|      0.0084432|    67|
| occupational therapy              |  0.0023249|     0.0080541|      0.0081471|    68|
| informatics                       |  0.0099826|     0.0080403|      0.0078628|    69|
| electrical engineering            |  0.0035042|     0.0078939|      0.0078633|    70|
| automation                        |  0.0022552|     0.0078823|      0.0082833|    71|
| geotechnology                     |  0.0027123|     0.0078614|      0.0080759|    72|
| engineering design                |  0.0024121|     0.0078412|      0.0081312|    73|
| forensics                         |  0.0026537|     0.0078035|      0.0076939|    74|
| education                         |  0.0047743|     0.0077080|      0.0075149|    75|
| mycology                          |  0.0025491|     0.0076594|      0.0076372|    76|
| logic                             |  0.0024510|     0.0076156|      0.0078584|    77|
| bibliometrics                     |  0.0027331|     0.0075872|      0.0077041|    78|
| anatomy                           |  0.0024009|     0.0075638|      0.0076517|    79|
| information systems               |  0.0036066|     0.0075381|      0.0072705|    80|
| electronics                       |  0.0037840|     0.0075347|      0.0074561|    81|
| Brazilian healthcare              |  0.0024541|     0.0075112|      0.0073808|    82|
| social work                       |  0.0028226|     0.0074891|      0.0074850|    83|
| biomaterials                      |  0.0112502|     0.0073331|      0.0067864|    84|
| logistics                         |  0.0034656|     0.0072290|      0.0071273|    85|
| human geography                   |  0.0045769|     0.0072145|      0.0071259|    86|
| ecology                           |  0.0061912|     0.0069880|      0.0068437|    87|
| astrophysics                      |  0.0042515|     0.0069680|      0.0064905|    88|
| nuclear engineering               |  0.0024137|     0.0069329|      0.0073273|    89|
| radioactivity                     |  0.0031276|     0.0069110|      0.0070901|    90|
| law reviews                       |  0.0026354|     0.0068964|      0.0067880|    91|
| mathematics                       |  0.0086304|     0.0068269|      0.0069352|    92|
| learning disabilities             |  0.0031019|     0.0067851|      0.0066117|    93|
| marine engineering                |  0.0023518|     0.0066475|      0.0070585|    94|
| geology                           |  0.0139842|     0.0064870|      0.0061994|    95|
| social justice                    |  0.0024211|     0.0063083|      0.0060054|    96|
| veterinary medicine               |  0.0057741|     0.0062613|      0.0060098|    97|
| spectroscopy                      |  0.0107509|     0.0061040|      0.0058136|    98|
| tourism                           |  0.0025631|     0.0060081|      0.0060159|    99|
| tribology                         |  0.0037282|     0.0059347|      0.0058167|   100|
| measurement                       |  0.0025076|     0.0057733|      0.0061204|   101|
| engineering management            |  0.0025191|     0.0054657|      0.0055476|   102|
| environmental science             |  0.0197079|     0.0054299|      0.0051928|   103|
| food science                      |  0.0125544|     0.0054239|      0.0051509|   104|
| chemistry & physics               |  0.0441748|     0.0053294|      0.0049305|   105|
| robotics                          |  0.0044990|     0.0049737|      0.0050510|   106|
| material science                  |  0.0073490|     0.0049279|      0.0049302|   107|
| agronomy                          |  0.0024144|     0.0048210|      0.0050570|   108|
| energy                            |  0.0072781|     0.0047346|      0.0046635|   109|
| textiles                          |  0.0023452|     0.0047212|      0.0052096|   110|
| wood                              |  0.0024436|     0.0045605|      0.0049897|   111|
| particle physics                  |  0.0051697|     0.0039993|      0.0039264|   112|

``` r
lastplace <- max(field_ranks$rank)
interestingfields <- c('mathematics',
                       'informatics',
                       'biomedical sciences',
                       'medicine',
                       'chemistry & physics',
                       'zoology',
                       'textiles')

library(ggrepel)
ggplot(field_ranks) +
  aes(rank, 100*BradleyTerry, label = field) +
  geom_point(colour = 'tomato2', size = 1) +
  geom_text_repel(data = subset(field_ranks, rank < 5 | field %in% interestingfields),
                  nudge_y = .15,
                  nudge_x = -5,
                  segment.alpha = .25,
                  family = 'Gill Sans MT Condensed', 
                  fontface = 'bold',
                  colour = 'tomato2',
                  point.padding = unit(0.1, 'lines')
                  ) +
  scale_x_reverse(name = NULL,
                  labels = scales::ordinal,
                  breaks = c(1, seq(20, lastplace, by = 20)),
                  minor_breaks = c(1, seq(5, lastplace, by = 5)),
                  limits = c(lastplace, 1),
                  expand = c(0.02, 0)) +
  scale_y_continuous(NULL,
                     position = 'right') +
  theme_bw() +
  theme(text = element_text(family = 'Gill Sans MT Condensed', face = 'bold'),
        axis.text = element_text(colour = 'tomato2'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = 'tomato2'),
        axis.line = element_line(colour = 'tomato2'),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.background = element_rect(fill = 'transparent', colour = NA))
```

``` r
interestingjournals <- c('R Journal',
                         'Stat Science',
                         'Biometrika',
                         #'J Stat Software',
                         'Annals',
                         'Stats in Medicine',
                         #'Biostatistics',
                         'Machine Learning',
                         'JRSS-B',
                         'Statistics',
                         #'J Applied Stats',
                         'Metrika',
                         'JCGS')
lastplace2 <- max(stats_ranks$rank)

ggplot(stats_ranks) +
  aes(rank, 100*BradleyTerry, label = journal) +
  geom_point(colour = '#2D7B95', size = 1) +
  geom_text_repel(data = subset(stats_ranks, rank <= 3 | journal %in% interestingjournals),
                  nudge_y = .15,
                  nudge_x = -5,
                  segment.alpha = .25,
                  family = 'Gill Sans MT Condensed',
                  fontface = 'bold',
                  colour = '#2D7B95',
                  point.padding = unit(0.2, 'lines')
                  ) +
  scale_x_reverse(name = NULL,
                  labels = scales::ordinal,
                  breaks = c(1, seq(20, lastplace2, by = 20)),
                  minor_breaks = c(1, seq(5, lastplace2, by = 5)),
                  limits = c(lastplace2, 1),
                  expand = c(0.02, 0)) +
  scale_y_continuous(NULL,
                     position = 'right') +
  theme_bw() +
  theme(text = element_text(family = 'Gill Sans MT Condensed', face = 'bold'),
        axis.text = element_text(colour = '#2D7B95'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = '#2D7B95'),
        axis.line = element_line(colour = '#2D7B95'),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.background = element_rect(fill = 'transparent', colour = NA))
```
