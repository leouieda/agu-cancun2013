# AGU Meeting of the Americas 2013

Code, text, notebooks, slides, etc for the AGU 2013 Meeting of the Americas.

# Abstract

## 3D magnetic inversion by planting anomalous densities

3D geophysical inversion
is usually computationally intensive.

Specially so are methods
that use a grid of rectangular prisms
for parametrization.

We present 
a 3D magnetic inversion algorithm
based on 
the computationally efficient method
of planting anomalous densities.

The algorithm consists
of an iterative growth
of the anomalous bodies
around prismatic elements 
called "seeds".

These seeds are user-specified
and have known magnetizations.

Thus, 
the seeds provide 
a way for the interpreter
to specify the desired skeleton
of the anomalous bodies.

The inversion algorithm
is computationally efficient 
due to various optimizations
made possible by 
the iterative nature
of the growth process.

The control provided 
by the use of seeds
allows one 
to test different hypothesis
about the geometry
and magnetization
of targeted anomalous bodies.

To demostrate this capability,
we applied our inversion method
to the Morro do Engenho (ME)
and A2 magnetic anomalies, Brazil 
(Figure 1a).

ME is an outcropping akaline intrusion
formed by dunites, peridotites and pyroxenites.

A2 is a magnetic anomaly 
to the Northeast of ME
and is thought to be
a similar intrusion.

ME has known magnetization,
therefore,
a plausible hypothesis 
is that A2 has the same magnetization.

We tested this hypothesis
by performing an inversion
using a single seed 
for each body.

Both seeds had 
the same magnetization.

Figure 1b shows that
the inversion produced
residuals up to 2000 nT
over A2
(i.e., a poor fit)
and less than 400 nT
over ME
(i.e., an acceptable fit).

Figure 1c shows that
ME is a compact outcropping body
with bottom at 
approximately 5 km.

However,
the estimate produced 
by the inversion
for A2 is outcropping
and is not compact.

Thus, 
the estimate for A2
provides a poor fit 
to the observations
and is not in accordance
with the geologic information.

This leads to the conclusion
that A2 does not have
the same magnetization as ME.
