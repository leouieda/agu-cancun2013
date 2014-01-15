Code, text, notebooks, slides, etc for the AGU 2013 Meeting of the Americas.

All **results** in the presentation were generated by the
[IPython notebooks](http://pyvideo.org/video/1744/teaching-with-the-ipython-notebook)
in the `notebooks` folder:

* [Results for the abstract](http://nbviewer.ipython.org/urls/raw.github.com/leouieda/agu-cancun2013/master/notebooks/abstract-results.ipynb)
* [Inversion of the ME and A2 anomalies](http://nbviewer.ipython.org/urls/raw.github.com/leouieda/agu-cancun2013/master/notebooks/ME-A2-inversion.ipynb)
* [Synthetic test with 2 bodies](http://nbviewer.ipython.org/urls/raw.github.com/leouieda/agu-cancun2013/master/notebooks/two-bodies.ipynb)
* [Run an inversion step-by-step](http://nbviewer.ipython.org/urls/raw.github.com/leouieda/agu-cancun2013/master/notebooks/steb-by-step.ipynb)
* [Extract seeds from a previous solution](http://nbviewer.ipython.org/urls/raw.github.com/leouieda/agu-cancun2013/master/notebooks/spine.ipynb)

Slides in PDF are available on figshare:
[doi:10.6084/m9.figshare.703651](http://dx.doi.org/10.6084/m9.figshare.703651)

Citation:

Uieda, L. and V. C. F. Barbosa (2013), 3D magnetic inversion by planting
anomalous densities, AGU Meeting of the Americas, Abstract GP22A-01.


# 3D magnetic inversion by planting anomalous densities

**Leonardo Uieda and Valéria C. F. Barbosa**

We present a new
3D magnetic inversion algorithm
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
To demonstrate this capability,
we applied our inversion method
to the Morro do Engenho (ME)
and A2 magnetic anomalies, central Brazil
(Figure 1a).
ME is an outcropping alkaline intrusion
formed by dunites, peridotites and pyroxenites
with known magnetization.
A2 is a magnetic anomaly
to the Northeast of ME
and is thought to be
a similar intrusion
that is not outcropping.
Therefore,
a plausible hypothesis
is that A2 has
the same magnetization
as ME.
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
approximately 5 km,
which is in agreement
with previous interpretations.
However,
the estimate produced
by the inversion
for A2 is outcropping
and is not compact.
In summary,
the estimate for A2
provides a poor fit
to the observations
and is not in accordance
with the geologic information.
This leads to the conclusion
that A2 does not have
the same magnetization as ME.
These results
indicate the usefulness
and capabilities
of the inversion method
here proposed.

**Figures**

![Figure1](https://raw.github.com/leouieda/agu-cancun2013/master/abstract-figs/me-a2.png)
Figure 1: a) total field magnetic anomaly map with the ME and A2 anomalies.
b) map of the inversion residuals (observed minus predicted data).
c) the two bodies estimated by the inversion.

**References**

Dutra, A. C., and Y. R. Marangoni (2009), Gravity and magnetic 3D inversion of
Morro do Engenho complex, Central Brazil, Journal of South American Earth
Sciences, 28(2), 193-203, doi:10.1016/j.jsames.2009.02.006.

Uieda, L., and V. C. F. Barbosa (2012a), Robust 3D gravity gradient inversion
by planting anomalous densities, Geophysics, 77(4), G55-G66,
doi:10.1190/geo2011-0388.1.
[[pdf](http://www.fatiando.org/papers/Uieda,Barbosa_2012(2).pdf)]


Uieda, L., and V. C. F. Barbosa (2012b), Use of the "shape-of-anomaly" data
misfit in 3D inversion by planting anomalous densities, SEG Technical Program
Expanded Abstracts, 1-6, doi:10.1190/segam2012-0383.1.
[[pdf](http://www.fatiando.org/papers/Uieda,Barbosa_2012(3).pdf)]
