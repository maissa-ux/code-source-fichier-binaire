
RDKit Enumeration Toolkit

RDKit Reaction Enumeration Toolkit tutorial.

Here you will learn how to enumerate reactions with various building blocks.

from __future__ import print_function
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.AllChem import ReactionFromRxnBlock, ReactionToRxnBlock
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

rxn_data = """$RXN

      ISIS     090220091539

  2  1
$MOL

  -ISIS-  09020915392D

  2  1  1  0  0  0  0  0  0  0999 V2000
   -2.0744    0.1939    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5440   -0.1592    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  1 F    2  17  35
V    1 halogen
M  RGP  1   2   1
M  ALS   1  2 F Cl  Br   
M  END
$MOL

  -ISIS-  09020915392D

  2  1  0  0  0  0  0  0  0  0999 V2000
    2.8375   -0.2500    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
    3.3463    0.0438    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
  1  2  1  0  0  0  0
V    2 amine.primary
M  RGP  1   1   2
M  END
$MOL

  -ISIS-  09020915392D

  3  2  0  0  0  0  0  0  0  0999 V2000
   13.5792    0.0292    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   14.0880    0.3229    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   13.0704    0.3229    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  RGP  2   2   1   3   2
M  END"""

rxn = ReactionFromRxnBlock(rxn_data)
rxn



Sanitizing Reaction Blocks

Reaction blocks come from many different sketchers, and some don't follow the MDL conventions very well. It is always a good idea to sanitize your reaction blocks first. This is also true for Smiles reactions if they are in kekule form.

AllChem.SanitizeRxn(rxn)


Preprocessing Reaction Blocks

You will note that there are some special annotations in the reaction block:

V   1 halogen
V   2 amine.primary

These allows us to specify functional groups with very specific smarts patterns.
These smarts patterns are preloaded into the RDKit, but require the use of PreprocessReactions to embed the patterns.

rxn.Initialize()
nWarn, nError, nReactants, nProducts, labels =  AllChem.PreprocessReaction(rxn)
print ("Number of warnings:", nWarn)
print ("Number of preprocessing errors:", nError)
print ("Number of reactants in reaction:", nReactants)
print ("Number of products in reaction:", nProducts)
print ("Preprocess labels added:", labels)




So now, this scaffold will only match the specified halogens and a primary amine. Let's get some!

!wget http://www.sigmaaldrich.com/content/dam/sigma-aldrich/docs/Aldrich/General_Information/1/sdf-benzylic-primary-amines.sdf -O amines.sdf

!wget http://www.sigmaaldrich.com/content/dam/sigma-aldrich/docs/Aldrich/General_Information/1/sdf-alkyl-halides.sdf -O halides.sdf


reagents = [
            [x for x in AllChem.SDMolSupplier("halides.sdf")],
            [x for x in AllChem.SDMolSupplier("amines.sdf")]
        ]

print ("number of reagents per template:", [len(x) for x in reagents])



Basic Usage
Creating a library for enumeration

Using the enumerator is simple, simply supply the desired reaction and reagents. The library filters away non-matching reagents by default. The RDKit will log any removed reagents to the info log.

library = rdChemReactions.EnumerateLibrary(rxn, reagents)


f you only want each reactant to match once ( and hence only produce one product per reactant set ) you can adjust the parameters:

params = rdChemReactions.EnumerationParams()
params.reagentMaxMatchCount = 1
library = rdChemReactions.EnumerateLibrary(rxn, reagents, params=params)



Enumerating the library

A library has an enumerator that determines what reagents are selected for purposes of enumeration. The default enumerator is a CartesianProduct enumerator, which is a fancy way of saying enumerate everything. You can get hold this enumerator by using the GetEnumerator method.

enumerator = library.GetEnumerator()
print (enumerator)

print ("Possible number of permutations:", enumerator.GetNumPermutations())


Understanding results of enumerations

Each enumeration result may contain multiple resulting molecules. Consider a reaction setup as follows:

A + B >> C + D

There may be multiple result molecules for a number of reasons:

    The reactant templates (A and B) match a reagent multiple times. Each match has to analyzed to form a new product. Hence, the result has to be a vector of products.

    There me be multiple product templates, i.e. C+D as shown above where C and D are two different result templates. These are output in a result as follows:

    result = enumerator.next()

    result == [ [results_from_product_template1], 
                [results_from_product_template2], ... ]

    result[0] == [results_from_product_template1]
    result[1] == [results_from_Product_template2]

Because there may be multiple product templates specified with potentially multiple matches, iterating through the results to get to the final molecules isa bit complicated and requires three loops. Here we use:

    result for the result of reacting one set of reagents
    productSet for the products for a given product template
    mol the actual product

In many reactions, this will result in a single molecule, but the datastructures have to handle the full set of results:

   for result in enumerator:
      for productSet in results:
          for mol in productSet:


count = 0
totalMols = 0
for results in library:
    for productSet in results:
        for mol in productSet:
            totalMols += 1
    count += 1
print("Number of result sets", count)
print("Number of result molecules", totalMols)



How does the enumerator work?

As mentioned, you can make a copy of the current enumeration scheme using the GetEnumerator method. Lets make a copy of this enumerator by copying it using copy.copy(..), this makes a copy so we don't change the state of the Library.




import copy
enumerator = copy.copy(library.GetEnumerator())
print(enumerator)
test_enumerator = copy.copy(enumerator)


Let's play with this enumerator.

First: let's understand what the position means (this is the same as library.GetPosition)

list(test_enumerator.GetPosition())



What this means is make the product from reagents[0][111] and reagents[1][130]

reagents[0][111]

reagents[1][130]




This also appears to be the last product. So lets' start over.

library = rdChemReactions.EnumerateLibrary(rxn, reagents, params=params)
test_enumerator = copy.copy(library.GetEnumerator())
list(test_enumerator.GetPosition())




We can Skip to the 100th result

test_enumerator.Skip(100)
pos = list(test_enumerator.GetPosition())
print(pos)

reagents[0][pos[0]]


reagents[0][pos[1]]


Let's advance by one here and see what happens. It's no surprise that for the CartesianProduct strategy the first index is increased by one.

pos = test_enumerator.next()
print(list(pos))



Enumeration States

Enumerations have states as well, so you can come back later using GetState and SetState

GetState returns a text string so you can save this pretty much anywhere you like.

Let's skip to the 100th sample and save both the state and the product at this step.


library = rdChemReactions.EnumerateLibrary(rxn, reagents, params=params)
# skip the first 100 molecules
library.GetEnumerator().Skip(100)
# get the state

state = library.GetState()
print("State is:\n", repr(state))

result = library.next()
for productSet in result:
    for mol in productSet:
        smiles = AllChem.MolToSmiles(mol)
        break




Now when we go back to this state, the next molecule should be the one we just saved.

library.SetState(state)
result = library.next()
for productSet in result:
    for mol in productSet:
        assert AllChem.MolToSmiles(mol) == smiles
        print(AllChem.MolToSmiles(mol), "==", smiles, "!")



Resetting the enumeration back to the beginning

To go back to the beginning, use Reset, for a CartesianProductStrategy this should revert back to [0,0] for indexing these reagents.

This is useful because the state of the library is saved when the library is serialized. See Pickling Libraries below.

library.ResetState()
print(list(library.GetPosition()))



Pickling Libraries

The whole library, including all reagents and the current enumeration state reagents is saved when the library is serialized.

s = library.Serialize() # XXX bug need default arg

library2 = rdChemReactions.EnumerateLibrary()
library2.InitFromString(s)

And the libraries are in lock step.

for i in range(10):
    result = library.next()
    for productSet in result:
        for mol in productSet:
            print("Result library1", AllChem.MolToSmiles(mol))
    result = library2.next()
    for productSet in result:
        for mol in productSet:
            print("Result library2", AllChem.MolToSmiles(mol))        




