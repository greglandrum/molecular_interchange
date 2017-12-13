# Header/disclaimers
**Draft date:** 13.12.2017

The ideas in here are not new. They are drawn from a bunch of conversations, both in person and online, that I've had over the past few years, my personal experience, and other versions of JSON for representing molecules. If you see something in here that looks like a wheel being reinvented, or that doesn't make sense to you, please let me know. I've picked JSON as a "physical representation" for the discussion just to make things concrete. The details of the actual format aren't super important, the flexible and extensible data model is the key part.

The immediate goal of this exercise is to end up with a useful chemical interchange format that can be used to facilitate transferring chemical information between the RDKit and other software. Ideally it will end up being adopted more broadly, but rather than doing the cat herding to start these conversations from a clean slate I would prefer to start from a working implementation.

# Design goals

## The problem we're attempting to solve

The primary goal here is to make it easier for software that's operating on molecules to interoperate. This need may come about because a "user" (I'll use that word loosely to cover users of scripting languages as well as developers of 3rd party software tools) wants to draw upon functionality from two different chemistry packages or because the developers of one package may want to use functionality from another package. The two ways that currently exist to do this are by using a molecule format like SDF/SMILES/etc or to directly write code that translates between the two different molecule representations. The first approach has problems because few (if any) of the existing formats have been designed with the goal of interoperability in mind. The second is painful because it likely requires directly linkage between the software packages and requires a separate translator for every toolkit pair.

## Requirements

### The format must be:
- efficient and easy to parse (this argues for JSON)
- versioned
- extensible in future versions (i.e. we should be able to add things that we miss in the initial design later)
- extensible by individual toolkits (see below)
- available in a text form
- open and well documented
- structured in such a way that toolkits can ignore pieces of information that they don't know how to/don't care to parse. Ideally this can be done without actually parsing that information
- it must have a simple and comprehensible chemistry model

### It would be nice if the format:
- supported a binary equivalent (messagepack for JSON, for example)
- were human readable/editable

### What should be representable
- discrete molecules, including biological macromolecules
- residue and chain information for macromolecules (analogous to info in the PDB)
- 0-N conformations per molecule
- partial atomic charges

### Out of scope for the core
- query information (at least for v1)
- molecules that aren't completely specified
- polymers/extended crystalline materials

## Extensibility by toolkits/software
In order to maximize efficiency, authors of toolkits may want to store precomputed data about a molecule into the interchange format. This could be used to store things like aromaticity, the results of ring perception algorithms, etc. Rather than try to come up with a set of these properties that make sense to everyone and then agree upon what they mean, we ensure that the interchange format is extensible in a way that allows toolkit authors to create blocks within the documents containing whatever information they want to pass along with the molecule.

# "Rules"
## Toolkit-specific blocks
- When parsing/writing an interchange document, toolkit-specific blocks that are not understood (or that haven't been parsed) should be passed along without modification.
- If the contents of the main block are changed at all, all other toolkit-specific blocks should be removed.
