#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

// std::endl to linebreak without /n
// std::vector<> er ish det samme som en liste i Python

//Kan lage selvkomponert header med funksjoner som du har deklarert og definert i andre .cpp. Dette kan du inkludere i andre koder senere.
//For å inkludere hpp-filer som du har laget selv i en kode må man compile med en path til filen. g++ -c main.cpp -I "navn på mappe til hpp"/
// g++ -c src/utils.cpp -I "navn på mappe til hpp"/
//For å lage .exe med hpp inkludert: g++ utils.o main.o -o main.exe
// Alt i ett: g++ main.cpp src/utils.cpp -I "navn på matte til hpp"/ -o main

//Type qualifier: "const int a = 1", const tells the code that the int variable a won't be modified.

//redirecting a terminal output: ./main.exe > output.dat. Hvis man kjører terminalen med det samme enda en gang, vil filen bare bli skrevet
//over med det samme.

//Using namespace std; Kan brukes for å slippe å skrive std::