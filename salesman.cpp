#include <cstdio>
#include <cstdlib>

#include <math.h>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <random>


using namespace std;

// simple structure to store city coordinates
// could also use std::pair<double> 
// or define a class

typedef struct {
  double lon, lat;
  //std::string name;
  char* name;
} COORD;

// functions expanded on below
// functions used to clean data
void setFilePath(char* fileName, char* input);
int getData(char* fname, COORD *cities);

// functions to execute the metropolis algorithm
void initialiseHot(std::vector<COORD>& route);
void sweep(std::vector<COORD>& route, double beta);
void reverseSection(std::vector<COORD>& route, int rand1, int rand2, int& rand3);
double calculateDistance(COORD cordinate1, COORD coordinate2);
double calculatePseudoEnergy(std::vector<COORD>& route);
void updateRoute();



int main(int argc, char *argv[]){

  // defines the variable in which the cities are stored
  const int NMAX = 2500;
  COORD cities[NMAX];

  // extracts information about cities from datafile
  char* inFileName = argv[1]; 
  setFilePath(inFileName, argv[1]); 
  std::cout << "inFileName:\t" << inFileName << std::endl;
  int nCities = getData(inFileName, cities);


  // vector of coordinates to be rearranged
  std::vector<COORD> route;
  std::vector<COORD> globalShortestRoute;
  //printf("Read %d cities from data file\n",ncity);
  //printf("Longitude  Latitude\n");
  for (int index = 0; index < nCities; index++){

    //printf("%lf %lf %s\n", cities[index].lon, cities[index].lat, cities[index].name);
    route.push_back(cities[index]);
  }




  // *craete output file with name: "salesman_cities.dat" by using outFileName = "salesman_" + inFileName
  std::string salesman_ = "salesman_";
  char* outFileName = strcat(salesman_.data(), inFileName);
  FILE *outFile = fopen(outFileName, "w");
  
  // begin the metropolis algorithm
  initialiseHot(route);
  double masterDistance = calculateDistance(route);

  unsigned nThermalisation = 100;
  unsigned nSweep = 50;

  int maxTemperature = 10;
  int nTemperature = 50;

  double meanPseudoEnergy;
  

  // loops through the prescribed number of temperatures
  for (unsigned indexTemperature = nTemperature; indexTemperature > 0; indexTemperature--) {
    
    initialiseHot(route); // rescrambles the route for every temperature

    // prescribes the correct temperature to this temperature step
    std::vector<COORD> tempShortestRoute = NULL; // stores the shortest route associated with a given temperature
    double temperature = (maxTemperature * indexTemperature) / nTemperature;
    if (temperature < 1e-6) temperature = 1e-6;  // protect agains overflow
    double beta = 1 / temperature;

    // reassigns the order of cities to let it adjust to the new temperature
    for (unsigned sweepIndex = 0; sweepIndex < nThermalisation; sweepIndex++){ sweep(route, beta); }

    meanPseudoEnergy = 0;
    
    // reassigns the order of cities for the purpose of minimising the distance
    for (unsigned sweepIndex = 0; sweepIndex < nSweep; sweepIndex++){

      sweep(route, beta);
      for (unsigned cityIndex = 0; cityIndex < (unsigned)nCities; cityIndex++){

	meanPseudoEnergy += 0;
	// *here we want to calculate the total path length of the route*
      }
      
      
    }

    // 
    tempShortestRoute = route;
    globalShortestRoute == NULL ? globalShortestRoute = tempShortestRoute;
    meanPseudoEnergy <= globalPseudoEnergy ? globalShortestRoute = tempShortestRoute;

    //fprintf();
  }

  //outFile->close();
  return 0;
}


// function sets the default path if none is found
void setFilePath(char* fileName, char* input){
  
  std::string defaultPath = "cities.dat";
  //input != NULL ? fileName = input: fileName = defaultPath.data(); std::cout << "function:\t" << fileName << std::endl;
  std::cout << "function:\t" << defaultPath.data() << std::endl;
}

// function reads city data from file 
int getData(char* fname, COORD *cities){
  
  std::cout << "function2:\t" << std::endl;
  //
  FILE* fp = fopen(fname,"r");
  const int bufsiz = 1000;
  char line[bufsiz + 1];
  int ncity = 0;
  
  //
  while(1){
    
    //
    fgets(line,bufsiz,fp);
    if (line[0] == '#') continue;  // skip comments
    if (feof(fp)) break;
    // we only scan for two numbers at start of each line
    sscanf(line,"%lf %lf %s", &cities[ncity].lon, &cities[ncity].lat, &cities[ncity].name);    
    ncity++;
  }
  fclose(fp);
  return ncity;
}

// calculates the distance between two points on the earth using the great circle method
double calculateDistance(COORD coordinate1, COORD coordinate2){

  double rE = 6371 * pow(10, 3);
  double deltaLat = coordinate2.lat - coordinate1.lat;
  double deltaLon = coordinate2.lon - coordinate1.lon;

  double a = pow(std::sin(0.5 * deltaLat), 2) + ( std::cos(coordinate1.lat) * std::cos(coordinate2.lat) * pow(std::sin(deltaLon), 2) );
  double c = 2*std::atan2(std::sqrt(a), std::sqrt(1 - a));

  return rE * c;
}


double calculatePseudoEnergy(std::vector<COORD>& route){

  double pseudoEnergy = 0;
  for (int index = 0; index < route.size() - 1; index++){

    pseudoEnergy += calculateDistance(route[index], route[index + 1]);
  }
  pseudoEnergy += calculateDistance(route[0], route[route.size() - 1]);

  return pseudoEnergy;
}

// randomly reorders the vector of cities to start
void initialiseHot(std::vector<COORD>& route){
  
  int size = route.size();
  for (int index = 0; index < size; index++){

    int temp = index + rand() % (size - index);
    std::swap(route[index], route[temp]);
  }
}

//
void sweep(std::vector<COORD>& route, double beta, double& masterDistance){
  
  int rand1, rand2;
  rand1 = rand() % route.size();
  do {
    
    rand2 = rand() % route.size();
  } while (rand1 == rand2 && (route.size() - rand1 + rand2) == 0);
  
  
  int rand3;
  vector<COORD> newRoute = reverseSection(route, rand1, rand2, rand3);
  double newDistance = masterDistance;

  double newDistanceSubtract = 0;
  double newDistanceAdd = 0;
  if (rand1 = 0 != true){

    newDistanceSubtract += calculateDistance(route[rand1], route[rand1 - 1]);
  } else {

    newDistanceSubtract += calculateDistance(route[0], route[route.size()]);
  }

  //
  if (rand2 = 0 != true){

    newDistanceSubtract += calculateDistance(route[rand2], route[rand2 - 1]);
  } else {

    newDistanceSubtract += calculateDistance(route[0], route[route.size()]);
  }

  //
  int groupSize = abs(rand1 - rand2);
  rand3 + groupSize > route.size() ? groupSize = route.size() - groupSize: ;
  if (rand3 = 0 != true){

    newDistanceAdd += calculateDistance(route[rand3], route[rand3 - 1]);
    newDistanceAdd += calculateDistance(route[rand3 + groupSize], route[rand3 + groupSize - 1]);

    rand1 < rand2 ? newDistanceAdd += calculateDistance(route[rand1 + groupSize], route[rand1 + groupSize - 1) : calculateDistance(route[rand2 + groupSize], route[rand2 + groupSize - 1) ;
  } else {

    newDistanceAdd += calculateDistance(route[0], route[route.size()]);
    newDistanceAdd += calculateDistance(route[groupSize], route[groupSize - 1]);
  }

  
  
  
  
  if ( newDistance <= masterDistance || drand48() - exp(beta * (newDistance - masterDistance))){

    masterDistance = newDistance;
    route = std::copy(newRoute.begin(), newRoute.end(), back_inserter(newRoute));
  }
}


vector<COORD> reverseSection(std::vector<COORD>& route, int rand1, int rand2, int& rand3){


  std::vector<COORD> routeCopy;
  std::copy(route.begin(), route.end(), back_inserter(routeCopy));
  if (rand1 - rand2 > 0) {

    rand3 = rand() % (n - abs(rand1 - rand2));
    std::vector<COORD> routeCopySection;
    std::copy(routeCopy[rand2], routeCopy[rand1], back_inserter(routeCopySection));
    routeCopy.std::erase(rand2, rand1);
    std::reverse(routeCopySection.begin(), routeCopySection.end());
    routeCopy.insert(routeCopy.begin() + rand3, routeCopySection.begin(), routeCopySection.end());
  } else {

    rand3 = rand() & abs(rand1 - rand2);
    std::vector<COORD> routeCopySectionLower;
    std::copy(routeCopy[0], routeCopy[rand1], back_inserter(routeCopySectionLower));

    std::vector<COORD> routeCopySectionHigher;
    std::copy(routeCopy[rand2], routeCopy[routeCopy.size() - 1], back_inserter(routeCopyHigher));

    routeCopy.std::erase(rand2, routeCopy.size() - 1); 
    routeCopy.std::erase(0, rand1);

    for (int index = 0; index < routeCopySectionLower.size(); index++){

      routeCopySectionHigher.push_back(routeCopySectionLower[index]);
    }

    std::reverse(routeCopySectionHigher.begin(), routeCopySectionHigher.end());
    routeCopy.insert(routeCopy.begin() + rand3, routeCopySectionHigher.begin(), routeCopySectionHigher.end());
  }

  return routeCopy;
}
