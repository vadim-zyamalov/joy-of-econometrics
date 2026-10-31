#import "@preview/fletcher:0.5.8" as fletcher: diagram, edge, node
#set text(size: 14pt)
#set page(height: auto, width: auto, margin: 1em)

#let nodes = ([$y_t$], [$y_(t-1)$], [$y_(t-2)$])
#let edges = ((1, 0), (2, 1), (2, 0))

#diagram(
  for (i, n) in nodes.enumerate() {
    node((i, 0), nodes.at(i), name: "t" + str(i))
  },
  for (from, to) in edges {
    let nf = label("t" + str(from))
    let nt = label("t" + str(to))
    edge(
      nf,
      nt,
      "->",
      bend: if calc.abs(from - to) > 1 { -30deg } else { 0deg },
      stroke: if calc.abs(from - to) > 1 { red } else { black },
    )
  },
)
