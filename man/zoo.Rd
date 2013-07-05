\name{zoo}
\docType{data}
\alias{zoo}
\title{zoo data}
\description{
The zoo data is a set from the UCI Machine Learning Repository (\url{http://archive.ics.uci.edu/ml/}).
This dataset contains 16 attributes, and 7 animal classes. The first column gives as descriptve name for each case. The next 16 columns each correspond to one feature. The last column is the classification information. The classification is a breakdown of which animals are in which of the 7 types.
}
\usage{
data(zoo)
}

\format{
\tabular{rl}{
[,1] \tab name: unique for each case\cr
[,2] \tab hair: Boolean\cr
[,3] \tab feathers: Boolean \cr
[,4] \tab eggs: Boolean \cr
[,5] \tab milk: Boolean \cr
[,6] \tab airborne: Boolean \cr
[,7] \tab aquatic: Boolean\cr
[,8] \tab predator: Boolean\cr
[,9] \tab toothed: Boolean\cr
[,10] \tab backbone: Boolean\cr
[,11] \tab breathes: Boolean\cr
[,12] \tab venomous: Boolean\cr
[,13] \tab fins: Boolean\cr
[,14] \tab legs. Set of values: [0,2,4,5,6,8]\cr
[,15] \tab tail: Boolean\cr
[,16] \tab domestic: Boolean\cr
[,17] \tab catsize: Boolean\cr
[,18] \tab Class labels, integer values in range [1,7].
}}

\details{
There are 7 classes all together:
\enumerate{
\item aardvark, antelope, bear, boar, buffalo, calf, cavy, cheetah, deer, dolphin, elephant, fruitbat, giraffe, girl, goat, gorilla, hamster, hare, leopard, lion, lynx, mink, mole, mongoose, opossum, oryx, platypus, polecat, pony, porpoise, puma, pussycat, raccoon, reindeer, seal, sealion, squirrel, vampire, vole, wallaby,wolf 
\item chicken, crow, dove, duck, flamingo, gull, hawk, kiwi, lark, ostrich, parakeet, penguin, pheasant, rhea, skimmer, skua, sparrow, swan, vulture, wren 
\item pitviper, seasnake, slowworm, tortoise, tuatara 
\item bass, carp, catfish, chub, dogfish, haddock, herring, pike, piranha, seahorse, sole, stingray, tuna 
\item frog1, frog2, newt, toad 
\item flea, gnat, honeybee, housefly, ladybird, moth, termite, wasp 
\item clam, crab, crayfish, lobster, octopus, scorpion, seawasp, slug, starfish, worm
}
}
\source{
Forsyth, R. (1990). UCI Machine Learning Repository - Zoo Data Set [http://archive.ics.uci.edu/ml/datasets/zoo]. Irvine, CA: University of California, School of Information and Computer Science.
}