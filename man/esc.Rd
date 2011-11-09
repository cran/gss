\name{esc}
\alias{esc}
\title{Embryonic Stem Cell from Mouse}
\description{
    Data concerning mouse embryonic stem cell gene expression and
    transcription factor association strength.
}
\usage{data(esc)}
\format{
    A data frame containing 1027 genes with the following variables.
    \tabular{ll}{
        \code{y1} \tab Gene expression after 4 days.\cr
        \code{y2} \tab Gene expression after 8 days.\cr
	\code{y3} \tab Gene expression after 14 days.\cr
	\code{klf4} \tab Score of TFAS with KLF4.\cr
	\code{nanog} \tab Score of TFAS with NANOG.\cr
	\code{oct4} \tab Score of TFAS with OCT4.\cr
	\code{sox2} \tab Score of TFAS with SOX2.\cr
	\code{clusterID} \tab Cluster identification.\cr
    }
}
\references{
    Gu, C. and Ma, P. (2010), Nonparametric regression with
    cross-classified responses.  Manuscript.
    
    Cai, J., Xie, D., Fan, Z., Chipperfield, H., Marden, J., Wong,
    W. H., and Zhong, S. (2010), Modeling co-expression across species
    for complex traits: insights to the difference of human and mouse
    embryonic stem cells.  \emph{PLoS Computational Biology}, \bold{6},
    e1000707.

    Ouyang, Z., Zhou, Q., and Wong, W. H. (2009), chip-seq of
    transcription factors predicts absolute and differential gene
    expression in embryonic stem cells.  \emph{Proceedings of the
    National Academy of Sciences of USA}, \bold{106}, 21521--21526.
}
\keyword{datasets}
