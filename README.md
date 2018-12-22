# Biers: Bayesian Inference Engines for Rethinking Statistics #

> "Give me some beers and I'll give you an impromptu lecture about
> that some time."
> -- Richard McElreath


This repo contains the output of my efforts to learn Bayesian data
analysis as presented in the splendid texbook "Statistical Rethinking"
by Richard McElreath <http://xcelab.net/rm/statistical-rethinking/>.
In my opinion, this book (combined with Richard's YouTube lecture
videos) is a fantastic introduction to statistical data analysis.  It
blends history and philosophy of science with technical matter, in
order to train empirical scientists in the techniques, strengths,
weaknesses, and risks of statistical modeling.


The textbook teaches students how to use the R programming language to
implement Bayesian stats, but I decided to make my life substantially
more complicated.  I was curious how far I could get implementing the
systems from this text in the Racket dialect of Scheme. 

R has been around for some time, has focused specifically on data
analysis, and has acquired a broad, deep community and substantial
mindshare.  It's libraries are well-tuned to this work.  Racket is not
nearly as widely used for data analysis.  Luckily for me a few
trailblazers have gone this route in the past, leaving behind some
artifacts that I could use along the way.

A guide to the core files (in conceptual order):

  * urn.rkt - From Chapter 2: a first taste of Bayesian inference,
    wherein you reason about the contents of an urn by repeatedly
    drawing balls from it, recording their color, and replacing them.
	
  * urn-abuse.rkt - Trolling the urn model: what happens when you use
    the urn model to reason about idealized coin toss data (i.e. what
    happens when your model poorly matches up with your data)?
		
  * globe.rkt - From Chapter 2: moving to continuous parameter values,
	reason about the proportion of water covering the earth by tossing
	a globe into the air and recording whether your finger lands on
	land or water.
	
	Also Chapter 3: Use random samples drawn from your
	posterior (and prior) distribution and your generative model to
	interrogate the model after (and before) fitting it to data.
	
  * heights.rkt - From Chapter 4: fit adult human height data to a
    normal distribution model; [COMING SOON!] use a linear regression
    model to reason about how adult human weight can predict human
    height.
	

A guide to the helper files:

  * diy.rkt - Some Racket code that implements routines that are
    already available as libraries.  Factored out of the main code as
    I learned to read the manuals, but maybe useful for reference and
    learning. 
  
  * laplace-approx.rkt - Laplacian quadratic approximation of
    distributions using normal distributions.  In essence this is
    McElreath's map/quap function reimplemented in Racket.  Depends on
    nlopt Racket package (patched to fix some bugs, see pull
    requests), which in turn depends on the NLopt C library, as well
    as racket-ml.

	  * <https://github.com/jkominek/nlopt>
	  * <https://nlopt.readthedocs.io/en/latest/>
	  * <https://github.com/danking/racket-ml>

  * diff.rkt - Function differentiation tools.  This is primarily used
    to numerically approximate Hessians, which are used for Laplacian
    quadratic approximation.  Partly inspired by an example
    differentiation routine from the Racket documentation.
	
  * summaries.rkt - Functions to compute summary statistics, either
    from samples or from a grid-based posterior approximation

  * utils.rkt - Some utility functions that have arisen during
    developments.  Factored out to a separate file.
	
  * images/ - Some screenshots of using the code and some helpful plots.
      
Some random extra files:  not directly related to McElreath's text,
but here are other related probability and stats files that may be of
interest.

  * subjective-priors.rkt - Mostly inspired by The Theory that Would
    Not Die by Sharon Bertsch McGrayne.  From the appendix: in
    the absence of much data, the (subjective) prior dominates the
    data.  One source of skepticism about Bayesian data analysis.
	
  * dice.rkt - Inspired by the article Pascal and the Invention of
    Probability Theory by Oystein Ore
	* <https://www.jstor.org/stable/2309286>


For illustration purposes, some tweets related to this effort:

  * <https://twitter.com/rg9119/status/1048263967972384768>
  * <https://twitter.com/rg9119/status/1057372763445293057>
  * <https://twitter.com/rg9119/status/1066059133940465665>


# Dependencies #

Here is a unified list of dependencies.  I'm grateful to all of the
folks who have made their code available online.  I've applied some
fixes/patches to get things to work here, noted below.  I will
eventually submit pull requests to the relevant repos in the hopes of
making the code here easier to get running.

  * racket-ml - a set of tools by Dan King for machine learning.  I
    specifically use code here for a multivariate normal
    distribution.  I munged the code a bit to overcome bitrot (see
    pull requests at github).

    * <https://github.com/danking/racket-ml>
	  
  * Racket-nlopt - wraps the NLOpt C library for optimizing functions.
    Some patches made (see pull requests on github). Depends on the C
    library.  Available via `raco pkg install nlopt` or on github.
	
	* <https://github.com/jkominek/nlopt> 
    * <https://nlopt.readthedocs.io/en/latest/> 

  * Racket_NumericalMethods - I use the adaptive numerical integration
  routine in this package for some early examples of non-grid bayesian
  inference.  It did not scale very far, but was at the least
  educational, and a source of at least one quite precise easy-to-implement
  reference implementation to compare against grid approximation.
  Used successfully in globe.rkt.  Used unsuccessfully in heights.rkt
  (I didn't get back even one pdf value after many minutes of waiting,
  but that could be my error).

	* <https://github.com/mkierzenka/Racket_NumericalMethods>

  * For experimentation, I recommend Alex Knauth's debug language for
    tracing intermediate executions.  Available via 
	`raco pkg install debug` or on github.
	
	* <https://github.com/AlexKnauth/debug>
