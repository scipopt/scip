# the sassy preprocessor for symmetry detection
The sassy preprocessor is designed to shrink large, sparse graphs. Before giving a graph to an off-the-shelf symmetry detection solver (such as [bliss](https://users.aalto.fi/~tjunttil/bliss/index.html), [dejavu](https://www.mathematik.tu-darmstadt.de/dejavu), [nauty](https://pallini.di.uniroma1.it/), [saucy](http://vlsicad.eecs.umich.edu/BK/SAUCY/), [Traces](https://pallini.di.uniroma1.it/)), the graph is instead first handed to the preprocessor. The preprocessor shrinks the graph, in turn hopefully speeding up the subsequent solver.

Some technicalities apply, though: a hook for symmetries must be given to sassy (a `sassy_hook`), and symmetries of the reduced graph must be translated back to the original graph. The preprocessor can do the reverse translation automatically, by providing a special hook that is in turn given to the backend solver (see the examples below). The graph format used by the preprocessor is described below as well.

The preprocessor comes in the form of a C++ header-only library and uses some features of C++17. To achieve good performance the library should be compiled with an adequate optimization level enabled (we use `-O3` for benchmarks), as well as assertions disabled (i.e., by using the flag `NDEBUG`).

## The graphs
We provide an interface for the construction of vertex-colored graphs in the class `static_graph`. The graph must first be initialized (either using the respective constructor or using `initialize_graph`). For the initialization, the final number of vertices and edges must be given. The number of vertices and edges can not be changed. Then, using `add_vertex` and `add_edge`, the precise number of defined vertices and edges must be added. The `add_vertex(color, deg)` function requests a color and a degree. Both can not be changed later (unless the internal graph is changed manually). Note that the function always returns the numbers `0..n-1`, in order, as the indices of the vertices. The `add_edge(v1, v2)` function adds an undirected edge from v1 to v2. It is always required that v1 < v2 holds, to prevent the accidental addition of hyper-edges. An example creating a path of length 3 is given below.

	#include "sassy/preprocessor.h"

	...

	sassy::static_graph g;
	g.initialize_graph(3, 2); // 3 vertices, 2 edges
	const int v1 = g.add_vertex(0, 1);
	const int v2 = g.add_vertex(0, 2);
	const int v3 = g.add_vertex(0, 1);
	g.add_edge(v1, v2);
	g.add_edge(v2, v3);

The internal graph format of sassy is `sgraph`. It follows the format of nauty / Traces closely, which is described in great detail [here](https://pallini.di.uniroma1.it/Guide.html). The `sgraph` format does not include a coloring, instead, an integer array `col` is given in addition to the graph. The meaning is that vertex `v` is mapped to color `col[v]`.

## The hook
The hook is a function that the user provides to symmetry detection software, which is called whenever a symmetry is found. It is the way symmetries are returned to the user. To maximize performance, it is however crucial to implement the hook in a certain way, which we describe in the following.

The definition for `sassy_hook` is as follows:

	typedef const std::function<void(int, const int *, int, const int *)> sassy_hook;

Note that a hook has four parameters, `int n`, `const int* p`, `int nsupp`, `const int* supp`. The meaning is as follows. The integer `n` gives the size of the domain of the symmetry, or in simple terms, the number of vertices of the graph. The array `p` is an array of length `n`. The described symmetry maps `i` to `p[i]`.

Crucially, `nsupp` and `supp` tell us which `i`'s are interesting at all: whenever `p[i] = i`, we do not want to iterate over `i`. To enable this, the array `supp` tells us all the points where `p[i] != i`. In particular, `supp[j]` for `0 <= j < nsupp` gives us the j-th vertex where `p[supp[j]] != supp[j]`. Note that `nsupp` gives the size of `supp`.  In many applications, reading symmetries in this manner is crucial for adequate performance.

An example is provided below:

	void my_hook(int n, const int *p, int nsupp, const int *supp) {
		for(int j = 0; j < nsupp; ++j) {
			const int i = supp[j];
			// do something with p[i]
		}
	}

The function `my_hook` can be wrapped into a `std::function` object as follows:

	auto hook = sassy::sassy_hook(my_hook);

`hook` can then be used as shown in the examples below.

## Example using bliss

	#include "bliss/graph.hh"
	#include "sassy/preprocessor.h"
	#include "sassy/tools/bliss_converter.h"

	...

	sassy::static_graph g;

	// graph must be parsed into g here!

	// lets preprocess...
	sassy::preprocessor p;
	// hook is a sassy_hook callback function
	p.reduce(&g, &hook);

	// ...and then we give the graph to bliss: first, convert the graph
	bliss::Graph bliss_graph;
	convert_sassy_to_bliss(&g, &bliss_graph);

	// then call bliss
	bliss::Stats bliss_stat;
	bliss_graph.find_automorphisms(bliss_stat, sassy::preprocessor::bliss_hook, (void*) &p);

	// done!

Note that the `bliss_hook` uses the field `p.saved_hook` to call the user-defined `sassy_hook` (i.e., in the example above `hook`). This also holds for all other solvers described below.


## Example using nauty

	#include "sassy/preprocessor.h"
	#include "sassy/tools/nauty_converter.h"
	#include "nauty/naugroup.h"

	...

	sassy::static_graph g;

	// graph must be parsed into g here!

	// lets preprocess...
	sassy::preprocessor p;
	// hook is a sassy_hook callback function
	p.reduce(&g, &hook);

	// ...and then we give the graph to nauty: first, convert the graph
	sparsegraph nauty_graph;
	DYNALLSTAT(int, lab, lab_sz);
	DYNALLSTAT(int, ptn, ptn_sz);
	convert_sassy_to_nauty(&g, &nauty_graph, &lab, &lab_sz, &ptn, &ptn_sz);

	// then call nauty
	statsblk stats;
	DYNALLSTAT(int, orbits, orbits_sz);
	DYNALLOC1(int,  orbits, orbits_sz, nauty_graph.nv, "malloc");
	static DEFAULTOPTIONS_SPARSEGRAPH(options);
	options.schreier = true;
	options.defaultptn = false;
	options.userautomproc = sassy::preprocessor::nauty_hook;
	if(nauty_graph.nv > 0) {
		sparsenauty(&nauty_graph, lab, ptn, orbits, &options, &stats, NULL);
	}

	// clean up
	DYNFREE(lab, lab_sz);
	DYNFREE(ptn, ptn_sz);
	SG_FREE(nauty_graph);

	// done!

Note that the `nauty_hook` uses the static field `preprocessor::save_preprocessor` to access `p` again, which in turn accesses `p.saved_hook`. If multi-threading is used in this configuration, `preprocessor::save_preprocessor` should be changed to `thread_local`. This also holds for Traces.

The `convert_sassy_to_nauty` method allocates memory for the graph, `lab` and `ptn` using the respective macros of nauty. Freeing up the memory has to be handled by the user.

## Example using Traces

	#include "sassy/preprocessor.h"
	#include "sassy/tools/traces_converter.h"
	#include "nauty/traces.h"

	...

	sassy::static_graph g;

	// graph must be parsed into g here!

	// lets preprocess...
	sassy::preprocessor p;
	// hook is a sassy_hook callback function
	p.reduce(&g, &hook);

	// ...and then we give the graph to Traces: first, convert the graph
	sparsegraph traces_graph;
	DYNALLSTAT(int, lab, lab_sz);
	DYNALLSTAT(int, ptn, ptn_sz);
	convert_sassy_to_traces(&g, &traces_graph, &lab, &lab_sz, &ptn, &ptn_sz);

	// then call Traces
	statsblk stats;
	DYNALLSTAT(int, orbits, orbits_sz);
	DYNALLOC1(int,  orbits, orbits_sz, traces_graph.nv, "malloc");
	static DEFAULTOPTIONS_TRACES(options);
	options.schreier = true;
	options.defaultptn = false;
	options.userautomproc = sassy::preprocessor::traces_hook;
	if(nauty_graph.nv > 0) {
		Traces(&traces_graph, lab, ptn, orbits, &options, &stats, NULL);
	}

	// clean up
	DYNFREE(lab, lab_sz);
	DYNFREE(ptn, ptn_sz);
	SG_FREE(traces_graph);

	// done!

The `convert_sassy_to_traces` method allocates memory for the graph, `lab` and `ptn` using the respective macros of Traces. Freeing up the memory has to be handled by the user.

## Example using saucy

	#include "sassy/preprocessor.h"
	#include "sassy/tools/saucy_converter.h"
	#include "saucy/saucy.h"

	...

	sassy::static_graph g;

	// graph must be parsed into g here!

	// lets preprocess...
	sassy::preprocessor p;
	// hook is a sassy_hook callback function
	p.reduce(&g, &hook);

	// ...and then we give the graph to saucy: first, convert the graph
	saucy_graph _saucy_graph;
	int* colors = nullptr;
	convert_sassy_to_saucy(&g, &_saucy_graph, &colors);

	// then call saucy
	struct saucy_stats stats;
	if(g.v_size > 0) {
		struct saucy *s = saucy_alloc(_saucy_graph.n);
		saucy_search(s, &_saucy_graph, 0, colors, &sassy::preprocessor::saucy_hook, &p, &stats);
		saucy_free(s);
	}

	// clean up
	delete[] colors;
	delete[] _saucy_graph.edg;
	delete[] _saucy_graph.adj;

	// done!

I want to mention that I have also seen a saucy version that uses a slightly different graph format. In this version, `saucy_graph` contains another field `colors`. In order to translate to this format, we just need to additionally set `_saucy_graph.colors = colors`, and remove `colors` from the parameter list of `saucy_search`:

...
	saucy_search(s, &_saucy_graph, 0, &sassy::preprocessor::saucy_hook, &p, &stats);
	...

Again, the `convert_sassy_to_saucy` method allocates memory for the graph and `colors`, which has to be handled by the user.

## Work in progress
Note that this project is still being actively developed. I am happy to take suggestions, bug reports, ...
