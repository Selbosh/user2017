useR!2017 poster
================
By David Selby and David Firth, Department of Statistics, University of Warwick.

A poster for the useR!2017 conference in Brussels.

### Title

Ranking influential communities in networks

### Abstract

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

extrafont::loadfonts(device = 'win', quiet = TRUE)
ggraph(stats_layout) +
  geom_edge_fan0(alpha = .05, colour = 'tomato2') +
  geom_node_point(aes(size = PageRank), fill = 'tomato2', pch = 21, colour = 'white') +
  geom_node_text(aes(label = name), size = 3,
                 repel = TRUE,
                 family = 'Arial Narrow',
                 fontface = 'bold',
                 colour = 'tomato2',
                 segment.alpha = .2) +
  coord_fixed() +
  theme_graph() +
  theme(legend.position = 'none')
```

![](img/stats_viz_web-1.png)

And here is a ranking of these statistical journals, by Bradley--Terry score.

``` r
dplyr::arrange(
  data.frame(journal = V(stats_ig)$name,
             Scroogefactor = scrooge::Scroogefactor(stats_xtab),
             PageRank = V(stats_ig)$PageRank,
             BradleyTerry = scrooge::BTscores(stats_xtab)),
  desc(BradleyTerry)
)
```

| journal              |  Scroogefactor|   PageRank|  BradleyTerry|
|:---------------------|--------------:|----------:|-------------:|
| J R STAT SOC A STAT  |      0.0282376|  0.0103465|     0.0283302|
| AM STAT              |      0.0270228|  0.0057930|     0.0270249|
| J R STAT SOC B       |      0.0256999|  0.0454141|     0.0259910|
| CAN J STAT           |      0.0239206|  0.0083131|     0.0241060|
| ANN STAT             |      0.0236113|  0.1100324|     0.0240983|
| BIOMETRIKA           |      0.0228638|  0.0498885|     0.0230736|
| J R STAT SOC C-APPL  |      0.0210658|  0.0089601|     0.0213421|
| TECHNOMETRICS        |      0.0209600|  0.0139015|     0.0208967|
| STAT SCI             |      0.0209032|  0.0153774|     0.0206936|
| BIOSTATISTICS        |      0.0209463|  0.0158200|     0.0205330|
| LIFETIME DATA ANAL   |      0.0182751|  0.0058858|     0.0199515|
| J STAT SOFTW         |      0.0201867|  0.0162305|     0.0197136|
| J AM STAT ASSOC      |      0.0194846|  0.0835435|     0.0194837|
| STAT COMPUT          |      0.0199332|  0.0109201|     0.0193991|
| J COMPUT GRAPH STAT  |      0.0198372|  0.0190766|     0.0193987|
| SCAND J STAT         |      0.0188394|  0.0169175|     0.0182769|
| BIOMETRICS           |      0.0184227|  0.0371235|     0.0181814|
| J MACH LEARN RES     |      0.0178295|  0.0302504|     0.0172690|
| MACH LEARN           |      0.0164286|  0.0057000|     0.0171842|
| STAT METHODS MED RES |      0.0168654|  0.0072661|     0.0167750|
| BERNOULLI            |      0.0163526|  0.0165876|     0.0162539|
| INT STAT REV         |      0.0161866|  0.0050933|     0.0162027|
| STAT NEERL           |      0.0157788|  0.0038321|     0.0161538|
| J TIME SER ANAL      |      0.0164238|  0.0075043|     0.0160605|
| J STAT PLAN INFER    |      0.0166547|  0.0368382|     0.0159328|
| AUST NZ J STAT       |      0.0156540|  0.0046732|     0.0158607|
| STAT MODEL           |      0.0157907|  0.0047046|     0.0157993|
| ANN I STAT MATH      |      0.0157902|  0.0083345|     0.0155305|
| EXTREMES             |      0.0149530|  0.0057294|     0.0153994|
| BAYESIAN ANAL        |      0.0154541|  0.0114376|     0.0153939|
| TEST-SPAIN           |      0.0162304|  0.0075401|     0.0152143|
| J QUAL TECHNOL       |      0.0133428|  0.0078417|     0.0144551|
| PROBAB ENG INFORM SC |      0.0129198|  0.0038546|     0.0141473|
| BIOMETRICAL J        |      0.0149472|  0.0107017|     0.0141281|
| STAT PROBABIL LETT   |      0.0146409|  0.0227675|     0.0137136|
| STAT MED             |      0.0145437|  0.0454388|     0.0136752|
| J NONPARAMETR STAT   |      0.0137078|  0.0080410|     0.0135831|
| COMPUT STAT DATA AN  |      0.0146516|  0.0420382|     0.0134674|
| J MULTIVARIATE ANAL  |      0.0142435|  0.0318751|     0.0132759|
| ENVIRON ECOL STAT    |      0.0129528|  0.0046764|     0.0131239|
| ENVIRONMETRICS       |      0.0120375|  0.0072662|     0.0129543|
| SCAND ACTUAR J       |      0.0119172|  0.0049733|     0.0124851|
| STOCH MODELS         |      0.0112702|  0.0041752|     0.0123314|
| ANN APPL STAT        |      0.0135044|  0.0192067|     0.0120295|
| METHODOL COMPUT APPL |      0.0108873|  0.0045467|     0.0115274|
| METRIKA              |      0.0094645|  0.0055462|     0.0102697|
| STAT METHODOL        |      0.0088005|  0.0046021|     0.0092631|
| ESAIM-PROBAB STAT    |      0.0092176|  0.0036076|     0.0092150|
| COMMUN STAT-SIMUL C  |      0.0095515|  0.0060782|     0.0091717|
| J STAT COMPUT SIM    |      0.0092627|  0.0082753|     0.0091212|
| COMMUN STAT-THEOR M  |      0.0094081|  0.0118277|     0.0090699|
| SEQUENTIAL ANAL      |      0.0088882|  0.0036688|     0.0089746|
| ELECTRON J STAT      |      0.0097919|  0.0154410|     0.0089544|
| J AGR BIOL ENVIR ST  |      0.0089445|  0.0045417|     0.0089017|
| REVSTAT-STAT J       |      0.0086771|  0.0028658|     0.0083682|
| SORT-STAT OPER RES T |      0.0079263|  0.0021687|     0.0083124|
| STATISTICS           |      0.0079254|  0.0050826|     0.0082635|
| STAT INTERFACE       |      0.0087419|  0.0040690|     0.0079890|
| J APPL STAT          |      0.0082399|  0.0066011|     0.0078485|
| APPL STOCH MODEL BUS |      0.0067881|  0.0035229|     0.0077174|
| STAT PAP             |      0.0065067|  0.0060219|     0.0069310|
| QUAL ENG             |      0.0063973|  0.0043052|     0.0068593|
| PAK J STAT           |      0.0058241|  0.0040568|     0.0064576|
| INT J BIOSTAT        |      0.0073174|  0.0051220|     0.0063488|
| BRAZ J PROBAB STAT   |      0.0056185|  0.0025517|     0.0059057|
| PHARM STAT           |      0.0055106|  0.0043995|     0.0058420|
| J BIOPHARM STAT      |      0.0059426|  0.0077599|     0.0057170|
| R J                  |      0.0063932|  0.0026105|     0.0054341|
| ASTA-ADV STAT ANAL   |      0.0051192|  0.0027805|     0.0052610|
| STAT METHOD APPL-GER |      0.0044261|  0.0028058|     0.0050218|
| QUAL RELIAB ENG INT  |      0.0044372|  0.0060060|     0.0046339|
| COMPUTATION STAT     |      0.0048448|  0.0035324|     0.0046084|
| J KOREAN STAT SOC    |      0.0045534|  0.0032268|     0.0045248|
| STAT BIOPHARM RES    |      0.0038185|  0.0033177|     0.0044349|
| REV COLOMB ESTAD     |      0.0041089|  0.0023813|     0.0043742|
| ADV DATA ANAL CLASSI |      0.0032452|  0.0025590|     0.0036667|
| QUAL TECHNOL QUANT M |      0.0021391|  0.0022255|     0.0021176|

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

extrafont::loadfonts(device = 'win', quiet = TRUE)
ggraph(other_layout) +
  geom_edge_fan0(alpha = .01, colour = 'tomato2') +
  geom_node_point(aes(size = PageRank), fill = 'tomato2', pch = 21, colour = 'white') +
  geom_node_text(aes(label = name), size = 3,
                 repel = TRUE,
                 family = 'Arial Narrow',
                 fontface = 'bold',
                 colour = 'tomato2',
                 segment.alpha = .2) +
  coord_fixed() +
  theme_graph() +
  theme(legend.position = 'none')
```

![](img/plot_stats_others-1.png)
