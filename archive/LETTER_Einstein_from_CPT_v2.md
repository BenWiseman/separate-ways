# Einstein's equations from CPT, with nothing quantised

B. H. Wiseman\
*Linnet Labs, Sydney, Australia (independent researcher; no external funding)*

## Abstract

Jacobson's derivation of the Einstein equation from the Clausius relation assumes three things it
does not supply: a temperature at every local Rindler horizon, an entropy proportional to that
horizon's area with the same coefficient everywhere, and such a horizon at every point in every
null direction. We show that taking CPT to be a symmetry of the universe, and not only of its
laws, supplies all three. The entropy law comes from the entanglement of the state the symmetry
forces and not from Bekenstein-Hawking, so nothing is borrowed back from the theory being derived.
The metric is a fixed background at every step and no gravitational degree of freedom is
quantised. Nineteen properties of general relativity follow against four inputs, among them the
sign of Newton's constant, four large spacetime dimensions, and $\Lambda$ as an integration
constant. Four observations would end the fold, two of them cosmological measurements now
being made.

## Main

General relativity is usually treated as a classical theory awaiting quantisation. The field
equations are read as the limit of something deeper, and the task is to find what. We show that
those equations follow instead from a symmetry, with the metric a fixed background at every step
and no gravitational degree of freedom promoted to an operator. On this reading there is nothing
about the field equations that a quantum theory of gravity still has to supply.

The symmetry is CPT, taken to hold of the universe rather than only of its laws. A law-level CPT
relates one process to another. A universe-level CPT says the spacetime is folded: two sheets, and
an antilinear involution $\Theta$ carrying each to the other. Boyle, Finn and Turok took this of
the bang and found a cosmology with no inflaton [1,2]. What follows here is what the same
assumption does at a horizon.

Write $\Theta=J\circ P_\perp$, with $J$ the antilinear modular conjugation of a wedge and
$P_\perp$ the antipodal map of the transverse directions. Two properties carry everything that
follows. First, $\Theta^2=1$, so the map is an involution and every field splits into the average
of its two readings and their difference. Those are the classical and quantum variables of the
Keldysh formalism, in use since the 1960s as a bookkeeping device and here given a reason to
exist: the classical world is where the two sheets agree. Second, $d\Theta=-\mathrm{Id}$. That is
the property which survives a change of frame, and it is the one the derivation uses.

The first property needs a moment, because it is what makes "nothing quantised" a position
rather than an evasion. The fold does not abolish the divide between classical and quantum
behaviour. It locates it, and then fixes its weight. For a single mode, with
$\langle\Phi_c^2\rangle=\frac12[W(0)+W_{\rm cross}]$ and
$\langle\Phi_q^2\rangle=2[W(0)-W_{\rm cross}]$, the cross-sheet correlator evaluated at the
half-period shift derived below gives, from the KMS condition alone,

$$\frac{\langle\Phi_q^2\rangle}{\langle\Phi_c^2\rangle}=4\tanh^2\frac{\beta\omega}{4}.$$

Below the horizon temperature the quantum half switches off as $(\beta\omega)^2/4$. That is a
classical limit reached in temperature and not in $\hbar$, and it says the classical world is a
horizon seen from below its own temperature. The two halves carry equal weight at
$\beta\omega=2\ln3$. Moving the shift off a half period moves that number, so it belongs to the
fold rather than to thermality.

The same construction says when the classical world switches on, and ties that moment to a
particle mass. Both are functionals of the same occupation, so the production integral that
fixes the dark-matter mass at the bang also fixes when the two branches stop interfering. With
the mass going as $I^{-2/5}$ and the crossing time as $M_1^{-1}I^{-2/3}$, eliminating the
production integral between them leaves

$$t_{\rm dec}\propto M_1^{2/3},$$

with $-1+\tfrac53=\tfrac23$ and nothing approximate in the step, so the exponent carries no free
parameter and is exact rather than fitted. At the adopted mass the crossing falls at
$1.417\times10^{-32}$ s, and at $2.249\times10^{-32}$ s if a Majorana pair is counted. A particle
mass and the moment the two time directions stop interfering are set by the same gravitational
event, neither chosen to suit the other. The normalisation is
state dependent and the power is not, which is why we quote the power.

**What the derivation needs.** Jacobson obtains the Einstein equation by applying
$\delta Q=T\,dS$ to a local Rindler horizon through every point in every null direction [3].
Because that step carries the whole of what follows, we re-derive it rather than cite it.
Integrating Raychaudhuri's equation near the bifurcation surface, where the expansion and shear
are first order so their squares are second, gives $\theta=-\lambda R_{kk}$. The heat flux carries
the same $\int\lambda\,d\lambda\,dA$, and the surface gravity cancels between the two sides, which
is what turns a statement about one accelerated observer into a field equation. What remains is
$2\pi T_{kk}=\eta R_{kk}$ for every null $k$. That a symmetric tensor annihilating every null
vector is a multiple of the metric is linear algebra and we check it rather than assert it: the
constraint map built from sixty random null vectors has a one-dimensional kernel, its ninth and
tenth singular values differing by $3\times10^{15}$, and that kernel is $g_{ab}$ to
$9\times10^{-16}$. Timelike vectors leave no kernel at all, so nullity is doing the work.

Four things go into that argument. The Clausius relation is the thermodynamic postulate it is made
of. The other three are a temperature at every such horizon, an entropy proportional to the
horizon area with a universal coefficient $\eta$, and horizons of that kind everywhere. In 1995
all three were assumed. Here they are consequences, and the rest of this Letter is the three derivations and what they assume.

**The temperature.** At a bifurcate Killing horizon the fold's map acts as a shift of Killing time
by half the periodicity, so the cross-sheet correlator is the direct one at $t-i\beta/2$. That is
a thermofield double at $\tanh r=e^{-\beta\omega/2}$. No temperature was put in and one came out.
What the fold supplies beyond thermality is that the shift is half a period and not a whole one,
which makes the map the square root of the thermal transformation and lets primitivity pick the
fundamental period.

**The horizons.** In any orthonormal frame the reflection of a Rindler wedge is
$\mathrm{diag}(-1,-1,+1,+1)$ and the transverse antipode is $\mathrm{diag}(+1,+1,-1,-1)$. Neither
factor is frame-independent. Their composition $-\mathrm{Id}$ is. A single map is therefore the
fold of every Rindler wedge through its own fixed point, at every boost and every orientation.
Composing the two factors in two hundred randomly boosted and rotated frames returns
$-\mathrm{Id}$ to $3\times10^{-14}$ while the reflections themselves move by order unity. This is
the step that makes the result a field equation rather than one equation of state per observer,
and it is the fold's, not an assumption about which horizons exist.

**The entropy law, and why this is not a loop.** This is the load-bearing step and the one on
which the whole argument can be dismissed, so we state the danger before the answer. If
$S\propto A$ were taken from Bekenstein-Hawking it would have come out of general relativity and
been fed back into a derivation of general relativity, and the result would be worth nothing.

It is not taken from there. Tracing out one wedge of the thermofield double the fold forces leaves
occupation $\sinh^2r$, and with $\tanh r=e^{-\beta\omega/2}$ that is exactly the Bose factor, so
the entanglement entropy of a wedge is the thermal entropy of its modes,
$s=(1+n)\ln(1+n)-n\ln n$. A horizon's modes carry a transverse momentum, the transverse directions
are translation invariant, and the mode count of a patch is therefore proportional to its area.
Every ingredient is flat-space field theory with nothing gravitational in it.

Proportionality is half of what Jacobson needs. The same $\eta$ must appear at every horizon, and
$\beta=2\pi/\kappa$ differs from one horizon to the next, so the occupation visibly carries the
surface gravity. It cancels. In the observer's own proper frequency the WKB mode count below
$\omega$ carries $\kappa$ in three places, and the Jacobian $d\omega=\kappa\,d\Omega$ cancels the
$1/\kappa$ the density carries. Computed with $\kappa$ kept throughout, $S/A$ does not move to one
part in $10^9$ while $\beta$ runs from 251 to 0.025. The reason is the theorem from which the
fold's own map is read: the wedge metric in boost coordinates contains no $\kappa$, and
Bisognano-Wichmann fixes the modular temperature at $2\pi$ in boost time at every wedge [4], so
$\beta\omega=2\pi\Omega$ identically. Let the modular temperature depend on $\kappa$ instead and
$S/A$ spreads by a factor of 21.6.

One condition attaches, and we found it by trying to break the result. $S/A$ depends on the cutoff as $1/\epsilon^2$, so a cutoff allowed to track the surface
gravity as $1/\kappa$ moves the answer by a factor of 52. The coefficient is the same at every
horizon provided the ultraviolet cutoff is one proper length and not one per horizon, which is
what a cutoff is. The mode count itself imposes no such condition: weighting the transverse
density by $e^{-k\epsilon}$ or by $1/(1+k^2\epsilon^2)$ moves $S/A$ from 2.137 to 1.629 and 1.931
while leaving the $\kappa$ independence at four parts in $10^{10}$, so the cancellation belongs to
the boost structure and not to one prescription.

What is claimed here is weaker than induced gravity and is exactly what the derivation needs.
Induced gravity computes $1/G$ from a given matter content and cutoff, and stands or falls on
whether the species count comes out. We claim nothing of that kind. The argument needs the
coefficient to be the same at every horizon, which is the paragraph above, and takes its value
from experiment, which is what $8\pi G=2\pi/\eta$ then reads as. Weaker is not empty: because
$\eta$ is fixed by a modular temperature and a transverse mode count, and neither is cosmological,
there is no dial in the fold that could make Newton's constant run.

**Horizons away from the fold's fixed points.** The fold's own fixed points are isolated, so at a
generic point no global fold acts. What the balance needs there is weaker than a global symmetry:
a map whose differential is $-\mathrm{Id}$, and one that is an isometry to the order the balance is
computed at. The geodesic symmetry $\exp_pv\mapsto\exp_p(-v)$ supplies both at every point of every
spacetime. In normal coordinates the metric's quadratic term is built from $R$ and its cubic term
from $\nabla R$, so the symmetry is exact through second order and fails at third. Cartan's theorem
says the failure is real wherever $\nabla R\ne0$, so what has to be settled is how heavy it is.
Measured on a surface with $\nabla K\ne0$, the residual $G_{\mu\nu}(v)-G_{\mu\nu}(-v)$ has exponent
$2.9994$ against 3 with a coefficient matching the closed form to five figures; on a surface whose
cubic part is harmonic, so that the curvature varies while $\nabla K$ vanishes at the point, the
exponent moves to 5.04. What the residual tracks is $\nabla R$ at the point and nothing else. Set
against the second-order term the balance uses, the ratio is 1.4435 times the patch size and
vanishes with it. The field equation is that limit.

Parity gives a second reason the balance cannot feel the failure. A map with differential
$-\mathrm{Id}$ multiplies a rank-$n$ tensor by $(-1)^n$, and the balance contains $R_{kk}$,
$T_{kk}$ and a transverse area element, every one of them even. An odd-rank quantity would break
and the balance holds none.

**What comes out.** Running the balance on the fold's own horizons gives

$$G_{ab}+\Lambda g_{ab}=8\pi G\left(T_{ab}^{\rm matter}+T_{ab}^{\rm img}[g,\Theta]\right),$$

with $\Lambda$ an integration constant of the derivation and not a prediction of it. That
$\Lambda$ enters this way is itself a result. The divergence step leaves $f+\eta R/2$ constant
across the whole spacetime, and nothing local fixes an integration constant; a boundary condition
does. The sign of Newton's constant comes with it, since $8\pi G=2\pi/\eta$ and $\eta$ is an
entropy per unit area, positive at every frequency. Gravity attracts here because entanglement
entropy is positive, and a priori that sign was free.

The second source is the one new term, and it is what makes this a different theory rather than
general relativity with a shifted coupling. A local curvature polynomial would renormalise $G$ and
do nothing else. The image stress is not local: it depends on the world function between a point
and its fold image, which is a bi-scalar. On the Einstein static universe, where homogeneity makes
every local invariant constant, it runs with the separation to the image, and the best local fit
leaves a residual of 139 per cent of itself while the same fitter reproduces a genuinely local
stress exactly. It carries no free parameter. Outside a horizon it is suppressed by mechanisms
rather than by a theorem: a point and its image are spacelike separated there, so a massive field's
image correlator carries $e^{-md}$, with $md$ above $10^{16}$ for an electron at a stellar horizon.
The largest value the term takes anywhere outside a horizon is $3\times10^{-34}$ of the dark
energy. It acts only inside a black hole, in a shell 2.1 microns thick at a solar mass, which is
the one region where general relativity already predicts its own breakdown.

**The ledger.** Figure 1 counts what the fold returns against what it is given. Nineteen
properties of general relativity come out and four go in. Where one property ends and the next
begins is a bookkeeping choice, and the table in Appendix A is there so the count can be taken apart; the lines themselves are produced by a script, and a checking
pass fails if the table and that script disagree.

Two of the four inputs are what writing a metric theory means, a Lorentzian signature and one metric carrying both the horizons and the
matter action. The other two are the measured constants $G$ and $\Lambda$. Those two could not
have come out, and the reason is dimensional rather than a matter of effort: in $\hbar=c=1$ every
structural input carries mass dimension zero, and $\eta=1/4G$ has dimension two. A symmetry of
dimension zero cannot produce a scale. What it produces is dimensionless numbers, and it produces
eight: $w_0=-1$ and $w_a=0$ exactly, a two-body line at exactly half the dark-matter mass, the
exponent $2/3$ tying a decoherence time to that mass, the half-period ratio that makes a horizon
thermal, the equal-weight frequency $2\ln3$, four large dimensions and no others, and the sign of
$G$.

Two entries usually placed on the input side are not there. The weak equivalence principle comes
out, being what the eikonal limit of a field on the shared metric does, checked against the
Christoffel system to $2.4\times10^{-12}$ with the mass cancelling out of the path. Matter
conservation comes out, from Noether's second theorem applied to a covariant matter action. What
the ledger does not touch: special relativity is input rather than output, and nothing here bears
on the gauge structure of matter. Those are not gaps in the derivation. They are what it is a
derivation from.

**What would break it.** The fold is more exposed than a reformulation has any right to
be, because the same fold that returns the field equations also fixes matter. Its state at the
bang sits half in each branch of an avoided crossing, which floors how many particles are made;
the observed abundance turns that floor into a ceiling, and the dark-matter fermion weighs at most
$491.6\pm2.0$ PeV, with its two-body neutrino line at $245.8\pm1.0$ PeV [5]. Exact stabilisation
floors the neutrino-mass sum at 58.8 meV and commits to constant dark energy, which shuts the
$w_0w_a$ escape that relaxes the present cosmological bound to 163 meV.

Four observations would end it, and none needs a threshold to be argued over: a securely
measured variation in Newton's constant at any level, since there is no dial to absorb it; a
robust detection of evolving dark energy, since an integration constant cannot evolve; a
cosmological neutrino-mass sum below 58.8 meV; and an identified decay above the two-body
endpoint. The first two test the structural half derived here, and the second two the matter
sector that comes with it.

**What this changes.** The field equations here follow from the fold and one constant, where they
had followed from a temperature, an entropy law and a constant. What they follow as is general
relativity plus one non-local term which carries no free parameter and which no observation ever
made could have reached. The cosmological constant is relocated rather than derived. It is a fact
about the fold's global state and not about the local physics, and the two places it could have
come from are both closed: the turning available in the late universe is short at every finite
radius, and at the bang an image stress independent of the scale factor would be of order $M_1^4$,
some 81.4 orders above the observed value. Naming it that way is more use than calling it
impossible, because it says where to look.

We do not claim that quantum gravity is impossible, or that the graviton has no commutator. The
free transverse-traceless algebra keeps its own, and the ordinary weak-field mediation calculation
remains available to anyone who wants it. Our claim is narrower and harder to work around: the
field equations do not wait on a quantum theory of gravity. They are what a universe folded by CPT
has to satisfy, and deriving them never requires the metric to be an operator.

## Figure

![](fig_ledger.png)

**FIG. 1.** The ledger. On the left is everything the fold is given; on the right is
everything it returns. Nineteen properties of general relativity come out against four inputs,
two of which are what writing a metric theory means and two of which are the measured constants.
Full wording of every line is in Appendix A.

## Acknowledgements

I thank Debbie Guimarães, Will Gittoes, Julie Wood, Emily Fountain, Sophie Bordson, and Jared
Young, and Nita Mannering who did not let me give up, and Angela Spina, who encouraged me to
follow this when the idea of past and future as folds first arrived, and who has backed it at
every point since. I thank Patricia Karr for her love, support, and always having my back. Ambre
Hammond for being the best cheerleader I could ask for. Joe Hanna for encouraging me to get back
into writing and publishing. I also thank my friends at Kandi Luxe, who have listened to more
cosmology over a cocktail bar than anyone signed up for: Jamila, Te, Gigi, Tuna, Layla, Claudia,
James, Jay, Michelle, Vesna, Amy, Adam, Mish, Bri, Hetti, Mehreen, Charlie, Petros, Blue, Noah,
Leah, Leyre, Hannah, and everyone else who has lovingly engaged with, or endured, my tangents and
rabbit holes. Aroha ahau ki a koutou katoa. And to Winston, fat and shameless as you are: there
is no more loyal configuration of matter in the cosmos. You are a good dog.

The author is responsible for every claim and decision here. Literature search, grammar,
proofreading, and automated tests, cross-checks, code revisions, and sanity checks used
contemporary research tools, large language models among them. No such tool is an author. Appendix C names
where the scripts behind these numbers live, and each check over them was validated by planting
an error it had to catch, since a check that matches nothing reports success.

## Competing interests

The author declares no competing interests. Linnet Labs is a company the author founded; it
neither funds nor directs this work, and there is no external funding.

## References

- [1] L. Boyle, K. Finn and N. Turok, Phys. Rev. Lett. **121**, 251301 (2018).
- [2] L. Boyle, K. Finn and N. Turok, Ann. Phys. **438**, 168767 (2022).
- [3] T. Jacobson, Phys. Rev. Lett. **75**, 1260 (1995).
- [4] J. J. Bisognano and E. H. Wichmann, J. Math. Phys. **17**, 303 (1976).
- [5] B. H. Wiseman, "Separate Ways and the Upside Down: the classical world and the field
  equations from a CPT fold, and the neutrino line that would break it", Zenodo (2026),
  doi:10.5281/zenodo.22888119.
- [6] C. Eling, R. Guedens and T. Jacobson, Phys. Rev. Lett. **96**, 121301 (2006).
- [7] G. Chirco and S. Liberati, Phys. Rev. D **81**, 024016 (2010).

## Appendices

### Appendix A. The ledger, line by line

Each derived line names the mechanism that produces it. Fuller wording, and the file behind each
number, are in Appendix E of ref. 5.

**The nineteen that come out.**

| No. | What comes out | From |
|---|---|---|
| 1 | The form of the field equations | Clausius on the fold's own horizons, Raychaudhuri integrated, the null-vector algebra checked, the Bianchi step written out |
| 2 | A temperature at every local Rindler horizon | The map is the half-period shift, so the cross-sheet correlator is a thermofield double |
| 3 | Those horizons, at every boost and orientation | The wedge reflection and the transverse antipode are each frame-dependent; their composition is not |
| 4 | An entropy proportional to horizon area | The entanglement entropy of the state the fold forces, with flat-space transverse mode counting. Bekenstein-Hawking is never invoked |
| 5 | The same entropy coefficient at every horizon, which makes the result Einstein's rather than one equation of state per observer | The surface gravity cancels between occupation and WKB mode density, because Bisognano-Wichmann fixes the modular temperature at $2\pi$ |
| 6 | The transverse involution at a black hole | The only free involutive isometry of a round bifurcation sphere is $-\mathrm{Id}$, and of Kerr's surface the antipodal map |
| 7 | That Newton's constant does not run | $\eta$ is fixed by a modular temperature and a transverse mode count, neither of them cosmological |
| 8 | The sign of Newton's constant | $8\pi G=2\pi/\eta$ and $\eta$ is an entropy density, which is positive |
| 9 | $\Lambda$ appearing at all, as an integration constant | The contracted Bianchi identity and matter conservation leave exactly one constant of integration |
| 10 | Null-focusing surfaces for Raychaudhuri to act on | Contact is conjugacy, so the fold's contact surfaces are caustics, at every charge and dimension |
| 11 | Four large spacetime dimensions and no others | The causal budget for contact closes only at $D=4$: the bill is $\pi$ and the interior pays $\pi/(D-3)$ per leg |
| 12 | The classical-quantum divide, and its amplitude | The fold's own parity, carrying weight $4\tanh^2(\beta\omega/4)$ |
| 13 | Black-hole thermodynamics | The same half-period shift, read at a bifurcate Killing horizon |
| 14 | The opposite time orientation of the two sheets | The modular flow of a wedge is the boost, positive on the right static patch and negative on its antipode |
| 15 | A hot bang, with no cold component at the fixed slice | $\Theta$ invariance needs $a$ odd in conformal time, which admits a fluid only where $-3w$ is an odd integer and excludes dust |
| 16 | The contracted Bianchi identity, which is what leaves exactly one integration constant | A theorem about the Levi-Civita connection of any metric; torsion breaks it at order unity |
| 17 | Matter conservation | Noether's second theorem on a covariant matter action, checked off shell in a generic curved metric |
| 18 | Matter following the metric's geodesics | The eikonal characteristics of a field on that metric are its geodesics, with the mass cancelling out of the path |
| 19 | Local Rindler horizons away from the fold's fixed locus | The geodesic symmetry has differential $-\mathrm{Id}$ and is an isometry through second order, which is the order the balance uses |

**The four that go in.**

| No. | What goes in | Why it could not come out |
|---|---|---|
| 1 | The value of Newton's constant | Every structural input has mass dimension zero and $\eta=1/4G$ has dimension two |
| 2 | The value of $\Lambda$ | An integration constant, and both places it could have come from are closed |
| 3 | Lorentzian signature | The fold keeps a parity exactly when the number of time directions is odd, which leaves $(1,3)$ and $(3,1)$ and excludes Euclidean |
| 4 | One metric, carrying both the horizons and the matter action | $T_{ab}$ means the variation of the matter action with respect to that metric and has no other meaning |

Two of the four are what writing a metric theory means and two are units.

### Appendix B. The two objections to the 1995 argument

Eling, Guedens and Jacobson showed the equilibrium Clausius relation fails once the entropy is not
proportional to area with a universal coefficient, and that an entropy-production term is needed
in its place; their worked case is $f(R)$ [6]. Here $\eta$ is not available to choose. It comes
out of the transverse mode count with the surface gravity cancelling, so the derivation sits in
the case the equilibrium relation was written for rather than assuming it does.

Chirco and Liberati showed the shear supplies an internal production term of its own, which they
identify with tidal heating [7]. That term is second order at the bifurcation surface, where the
expansion and the shear are both first order, and the Raychaudhuri measurement says how far out
that survives: the linear behaviour holds to $4\times10^{-4}$ at $\lambda=0.05$ and visibly worse
beyond. What remains assumed is the equilibrium reading itself, taken near the bifurcation surface
of each local wedge. That is the first of the four, and the one the fold was never going to
supply.

### Appendix C. Code and data

Every number quoted here is produced by a script in the release accompanying ref. 5, at
https://github.com/BenWiseman/separate-ways, each named for the claim it supports. A checking pass
fails if a quoted number and the script that produces it come apart.
