# Header/disclaimers
**Draft date:** 12.05.2018

The ideas in here are not new. They are drawn from a bunch of conversations, both in person and online, that I've had over the past few years, my personal experience, and other versions of JSON for representing molecules. If you see something in here that looks like a wheel being reinvented, or that doesn't make sense to you, please let me know. I've picked JSON as a "physical representation" for the discussion just to make things concrete. The details of the actual format aren't super important, the flexible and extensible data model is the key part.

The immediate goal of this exercise is to end up with a useful chemical interchange format that can be used to facilitate transferring chemical information between the RDKit and other software. Ideally it will end up being adopted more broadly, but rather than doing the cat herding to start these conversations from a clean slate I would prefer to start from a working implementation.

This document captures my thinking on what we're trying to accomplish. There's an actual spec that more or less aligns with this here: https://github.com/CommonChem/CommonChem

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
- multiple molecules
- discrete molecules, including biological macromolecules
- residue and chain information for macromolecules (analogous to info in the PDB)
- 0-N conformations per molecule
- partial atomic charges

### Out of scope for the core
- query information (at least for v1)
- molecules that aren't completely specified
- polymers/extended crystalline materials
- chemical reactions

## Extensibility by toolkits/software
In order to maximize efficiency, authors of toolkits may want to store precomputed data about a molecule into the interchange format. This could be used to store things like aromaticity, the results of ring perception algorithms, etc. Rather than try to come up with a set of these properties that make sense to everyone and then agree upon what they mean, we ensure that the interchange format is extensible in a way that allows toolkit authors to create blocks within the documents containing whatever information they want to pass along with the molecule.

# "Rules"
## Toolkit-specific blocks
- When parsing/writing an interchange document, toolkit-specific blocks that are not understood (or that haven't been parsed) should be passed along without modification.
- If the contents of the main block are changed at all, all other toolkit-specific blocks should be removed.

# Example(s)
These are expressed in JSON since it's a simple functional way to show the data model. As mentioned above, we aren't tied to only using JSON.

## Small molecule without coordinates
This corresponds to the SMILES `c1c(C=CC)cccc1O/C=C\\[C@H]([NH3+])Cl` and includes molecular properties:
```
{
  "commonchem": {
    "version": 10
  },
  "defaults": {
    "atom": {
      "z": 6, "impHs": 0, "chg": 0, "nRad": 0, "isotope": 0, "stereo": "unspecified"
    },
    "bond": {
      "bo": 1, "stereo": "unspecified"
    }
  },
  "molecules": [
    {
      "atoms": [
        { "impHs": 1 },
        {},
        { "impHs": 1 },
        { "impHs": 1 },
        { "impHs": 3 },
        { "impHs": 1 },
        { "impHs": 1 },
        { "impHs": 1 },
        {},
        { "z": 8 },
        { "impHs": 1 },
        { "impHs": 1 },
        { "impHs": 1, "stereo": "ccw" },
        { "z": 7, "impHs": 3, "chg": 1 },
        { "z": 17 }
      ],
      "bonds": [
        {
          "bo": 2,
          "atoms": [
            0, 1
          ]
        },
        {
          "atoms": [
            1, 2
          ]
        },
        {
          "bo": 2,
          "atoms": [
            2, 3
          ]
        },
        {
          "atoms": [
            3, 4
          ]
        },
        {
          "atoms": [
            1, 5
          ]
        },
        {
          "bo": 2,
          "atoms": [
            5, 6
          ]
        },
        {
          "atoms": [
            6, 7
          ]
        },
        {
          "bo": 2,
          "atoms": [
            7, 8
          ]
        },
        {
          "atoms": [
            8, 9
          ]
        },
        {
          "atoms": [
            9, 10
          ]
        },
        {
          "bo": 2,
          "atoms": [
            10, 11
          ],
          "stereo": "cis",
          "stereoAtoms": [
            9, 12
          ]
        },
        {
          "atoms": [
            11, 12
          ]
        },
        {
          "atoms": [
            12, 13
          ]
        },
        {
          "atoms": [
            12, 14
          ]
        },
        {
          "atoms": [
            8, 0
          ]
        }
      ],
      "extensions": [
        {
          "name": "rdkitRepresentation",
          "formatVersion": 1,
          "toolkitVersion": "2018.09.1dev1",
          "aromaticAtoms": [
            0, 1, 5, 6, 7, 8
          ],
          "aromaticBonds": [
            0, 4, 5, 6, 7, 14
          ],
          "cipRanks": [
            6, 8, 3, 1, 0, 4, 2, 5, 10, 13, 9, 7, 11, 12, 14
          ],
          "cipCodes": [
            [
              12, "R"
            ]
          ],
          "atomRings": [
            [
              0, 8, 7, 6, 5,1
            ]
          ]
        }
      ]
    }
  ]
}
```

## Small molecule with multiple conformers

Corresponds to SMILES `O[C@H]([35Cl])F` and includes partial charges:
```
{"commonchem": {"version": 10},
 "defaults": {
   "atom": {"z": 6, "impHs": 0, "chg": 0, "stereo": "unspecified", "nrad": 0},
   "bond": {"bo": 1, "stereo": "unspecified", "stereoAtoms": []}
 },
 "molecules": [{"name": "example 2", "atoms": [{"z": 8},
 {"z": 6, "stereo": "cw"},
 {"z": 17, "isotope": 35},
 {"z": 9},
 {"z": 1},
 {"z": 1}],
"bonds": [{"atoms": [0, 1]},
 {"atoms": [1, 2]},
 {"atoms": [1, 3]},
 {"atoms": [0, 4]},
 {"atoms": [1, 5]}],
"conformers": [{"dim": 2, "coords": [[-1.1988, -0.0452], [-0.3332, 0.4556], [0.166, 1.322], [0.5332, -0.0438], [-2.0652, 0.4542], [-0.834, 1.3212]]},
 {"dim": 3, "coords": [[0.9554, -0.3743, -0.4679], [-0.3616, -0.1311, -0.0795], [-0.5097, 1.6256, 0.1706], [-0.5779, -0.792, 1.1042], [1.5277, 0.1064, 0.1628], [-1.0339, -0.4346, -0.8902]]}],
"extensions": [{"name": "rdkitRepresentation", "toolkitVersion": "2018.03.1.dev1", "formatVersion": 1, "aromaticAtoms": [],
"aromaticBonds": [],
"cipCodes": [[1, "S"]],
"atomRings": []},
{"name": "partialCharges", "chargeType":"gasteiger", "formatVersion":10, "generator":"RDKit", "generatorVersion": "2018.03.1.dev1", "values": [-0.352, 0.273, -0.055, -0.198, 0.215, 0.117]}
]}]}
```


# Notes
- the usage of defaults objects (e.g. atomDefaults and bondDefaults) is an idea borrowed from the JSON format produce by the OpenEye tools. It allows files to be kept shorter (reducing the amount of parsing that has to be done) while still including all information in the file (i.e. you don't need the documentation to know that the default number of radicals on an atom is zero).
