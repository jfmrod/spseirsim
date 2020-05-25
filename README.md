SPSEIRSIM v1.0
Joao F Matias Rodrigues <jfmrod@gmail.com>

A fast multithreaded spatial, age and household structured epidemic SEIR model designed to use high resolution geographic estimates of demographic data.


# Description

The model uses as input age-structured demographics data from country census reports, as well as real country boundaries and high resolution estimated population counts from geographic databases.
Household structure is generated randomly following loosely reported household size distributions for the UK (also used in the model from Ferguson et al). Heuristic rules are used to generate realistic households.

Epidemiological parameters for SARS-Covid-2 were obtained from the model of Davis et al which are similar to the ones used by Ferguson et al.

The extreme age-dependence in mortality and hospitalization rates suggest that the herd immunity threshold can be reached in a population of 10 million in a few months with minimal fatalities and low hospitalization rates if only a younger segment of the population is exposed to the virus. Initial static analysis assuming the limit case of 100% exposed individuals of the different segments of the population shows that exposing only the sub 50 or 55 years old would yield a proportion of exposed population of 60% to 70%, large enough for herd immunity with a low number of hospitalizations and fatalities.
The model made available here was used to perform a preliminary dynamical analysis of this epidemiological scenario and to evaluate its feasibility. The questions addressed were: the impact of different levels of protection of the older higher risk population, identify policies needed to reach the herd immunity threshold (fraction of infected people upon which the removal of all isolation measures does not lead to a new epidemic wave) with the following constraints: 1) without exceeding hospital capacities and reducing overall fatalities, and 2) minimum isolation times and strength of isolation measures.


What is modelled:
- Individuals with specific age are assigned to each household for the whole population in a country (tested on up to 80 million people when simulating Germany, but should work for bigger countries)
- Households assigned to locations on the grid satisfying the constraints of the number of people living in that location
- Full state simulation: per household and per individual states
- Probability of transmission for each individual in a household is calculated based on:
  a) household: number of active infected (asymptomatic, presymptomatic, symptomatic) individuals in household
  b) local: position on a geographical grid (30 arc second) depending on the proportion of infected local and neighboring grid cells (using a gaussian kernel function)
  c) global: total proportion of infected individuals in the whole region
- Delays between states (exposed, asymptomatic, presymptomatic, symptomatic, time to hospitalization, time in ICU or NonICU, time to fatal outcome)
- Proportions between different states (asymptomatic/symptomatic, hospitalization rates, ICU/NonICU rates, fatality rates)
- Isolation measures reduce global and local transmission probabilities, but not intra-household transmission probabilities
- Age threshold dependent isolation measures

What is not modelled:
- Explicit international or interregional travel (rail, road, air), only a global factor is implemented
- Explicit travel, workplace, schools and other group activies, individuals become infected assuming they interact randomly with actively infected individuals in their local area
- Separate effect of isolation measures on strength and number of contacts and travel distances
- Regional isolation triggers
- Nursing homes and institutions housing larger numbers of people

The current implementation of the model would allow for the straightforward modelling of such effects, but a meaningful simulation incorporating such effects
would require more assumptions, further data and validation. To avoid having too much complexity therefore a conscious decision was made not to implement these effects explicitly.

Household age assignment:
Individual ages are assigned to households by matching two adults per household of similar ages (0 to 15 age difference).
Younger dependent children (up to 20 years old) are in a next step assigned to households of size 3 or larger with an age difference to adults of 25 to 35 years.



# Installation

Libraries:
- Required: libeutils (source available using svn from https://www.konceptfx.com/svn/eutils)
- Required: libshp libgeotiff libtiff
- Optional (for video generation): libavcodec libcairo

Required datasets (not included):
- World map population counts in Geotiff format (30 arc sec)
- Shape file for country
- Current age segmented population census (5 year segments)

World map population counts (GeoTiff) can be obtained from:
http://sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse

File used: gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif


Shape file for European countries in angular (lat/lon) coordinates
https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/countries

File used: NUTS_RG_01M_2016_4326_LEVL_2.shp

Current age structured population census:
https://github.com/cmmid/covid-uk.git

File used: covidm/data/wpp2019_pop2020.rds

Original data from: https://population.un.org/wpp/Download/Standard/Population/


Epidemiological parameters obtained from Davis et al. model (https://github.com/cmmid/covid-uk.git)

Distribution parameters:
- Incubation time                                  (Average: 4 days)   Gamma distributed (shape=4,scale=4)
- Asymptomatic time                                (Average: 5 days)   Gamma distributed (shape=5,scale=5)
- Presymptomatic time                              (Average: 1.5 days) Gamma distributed (shape=1.5,scale=4)
- Symptomatic time                                 (Average: 3.5 days) Gamma distributed (shape=3.5,scale=4)
- Time to hospitalization (NonICU) after symptoms: (Average: 7 days)   Gamma distributed (shape=7,scale=7)
- Time to ICU hospitalization after symptoms:      (Average: 7 days)   Gamma distributed (shape=7,scale=7)
- Time in hospital (NonICU)                        (Average: 8 days)   Gamma distributed (shape=8,scale=8)
- Time in ICU                                      (Average: 10 days)  Gamma distributed (shape=10,scale=10)
- Time to fatal outcome since symptomatic          (Average: 22 days)  Gamma distributed (shape=22,scale=22)

Proportions:
- Symptomatic/asymptomatic infections             0.66    2 symptomatic : 1 asymptomatic

Age dependent proportions: (can be found in data/agegroup.infparams)
- Symptomatic hospitalizations
- ICU/NonICU hospitalizations
- ICU fatalities
- NonICU fatalities 
 


2) Running a simulation:

# Parameters:
-fglobal        Global transmission factor, represents travel between simulated regions (not explicitly modelled)
-fintra         Intra-household transmission factor, affects spread of infections inside households
-flocal         Local transmission factor (radius determined by -ftravel), affects spread of infections in adjacent grid positions (in a 11x11 grid block)
-ftravel        Standard deviation of gaussian kernel for spread of infection in number of cells (dependent on cell units, should be angular units)
                width of kernel is corrected for average latitude of simulated region
-r0             Basic reproduction number
-iseed          Number of starting exposed individuals placed on grid with highest population count or custom lat:lon:count can be provided comma separated for multiple
                locations. Example: 20  or  46.0:8.95:10,46.2044:6.1432:10

-cf             Correction factor for IFR and hospitalization rates. New evidence based on serological studies indicates that previous estimates may be over estimated by almost 3 fold (see Streeck et al). Fatality rates for SARS-Covid-2 may however be underestimated by about 20% to 30% as country wide aggregated general mortality trends show mortality rates above the expected (compared to historical values) even when accounting for Covid19 deaths. These are probably fatal cases that have not reached the hospitals and were therefore not tested nor classified as covid-19 fatalities.


-tmax           Simulated time period (days)


-nthreads       Number of threads used for multithreading


Example: (No measures implemented)

./spseirsim -fglobal 0.1 -ftravel 4.0 -r0 2.7 -tmax 365 -iseed 20 -nthreads 4


# Video output: 

-ovideo      Specify output filename for video

Example:
./spseirsim -fglobal 0.1 -ftravel 4.0 -r0 2.7 -tmax 365 -nthreads 4 -iseed 20 -ovideo simulation.mp4


# With isolation measures triggered by new ICU counts:

-icutrigger        Number of new ICU cases per day that trigger isolation measures
-ficu              Initial isolation measure strength when triggered (0.0 perfect isolation, 1.0 no isolation)
-iculowlimit       ICU threshold at which isolation measures may be reduced (requires ICU cases to be reducing daily and -rthres threshold condition be met)
-rthres            Threshold of effective R value before isolation measures may be reduced
-fstep             Amount of reduction when isolation measures are reduced
-triggerdays       Number of days that must pass before new reduction in measures is considered
-triggercount      Number of reduction steps before all isolation measures are removed

Example:
./spseirsim -fglobal 0.1 -ftravel 4.0 -r0 2.7 -tmax 365 -nthreads 4 -iseed 20 -icutrigger 25 -rthres 0.95 -iculowlimit 1500 -ficu 0.25 --triggerdays 4 -fstep 0.10 --triggercount 5


# With age group specific isolation measure triggers:

-foiso             Strength of older group isolation measures
-oisotrigger       Trigger older group isolation on specified number of new cases (Symptomatic individuals)
-oisorelease       Trigger reduction of older group isolation on specified number of new cases (should be lower than -oisotrigger)
-oisorelease2      Trigger removal of older group isolation measures on number of new cases

Example:
./spseirsim -fintra 1.2 -finter 1.0 -fglobal 0.1 -ftravel 4.0 -r0 2.7 -tmax 365 -nthreads 40 -agethres 55 -iseed 46.0:8.95:10,46.2044:6.1432:10 -icutrigger 45 -rthres 1.0 -iculowlimit 1600 -ficu 0.5 --triggerdays 4 -fstep 0.1 -cf 0.37 -oisotrigger 4000 -oisorelease 6000 -oisorelease2 1000 -triggercount 1 -foiso 0.05


# With isolation measures activated on specified days:

-events        Sets isolation measure for all individuals (global and local reduction in transmission) at a certain days <day1>:<isolation_strength1>,<day2>:<isolation_strength2>,...
-oevents       Sets older population isolation measures on specific days  <day1>:<isolation_strength1>,<day2>:<isolation_strength2>,...

Example:
./spseirsim -fglobal 0.1 -ftravel 4.0 -r0 2.7 -tmax 365 -nthreads 4 -iseed 20 -events 10:0.3,120:0.4,150:1.0

This sets the isolation strength to 0.3 on day 10, then to 0.4 on day 120, and then to 1.0 (no isolation) on day 150.


# Specifying shape files, selecting the shapes or country(ies) when using the EU NUTS shape file:

-fshape        Specify shape file
-fpop          Specify age structured population counts
-shpsel        Specify selection rule for shapes in shape file, format:    FIELD1=VALUE1,FIELD2=VALUE2 or FIELD1=VALUE1,-FIELD2=VALUE2 (FIELD2=VALUE2 entries are excluded)

Example:
./spseirsim -tmax 365 -iseed 20 -fpop data/portugal.agegroups -shpsel 'CNTR_CODE=PT,-NUTS_ID=PT20,-NUTS_ID=PT30'

Selects regions belonging to continental Portugal but excludes non-continental regions (Madeira and Acores islands)
Fields in shape files can be inspected using the dbfdump utility from the shapefile tools.

Example of population age groups file format:

"country_code"	"name"	"age"	"f"	"m"	"location_type"
620	"Portugal"	"0-4"	194.288	206.364	4
620	"Portugal"	"5-9"	213.738	226.686	4
620	"Portugal"	"10-14"	239.813	250.186	4
620	"Portugal"	"15-19"	256.708	268.459	4
620	"Portugal"	"20-24"	269.888	269.946	4
620	"Portugal"	"25-29"	270.969	262.895	4
620	"Portugal"	"30-34"	287.271	273.726	4
620	"Portugal"	"35-39"	336.709	317.603	4
620	"Portugal"	"40-44"	395.961	369.49	4
620	"Portugal"	"45-49"	420.341	390.119	4
620	"Portugal"	"50-54"	389.162	354.227	4
620	"Portugal"	"55-59"	390.593	347.025	4
620	"Portugal"	"60-64"	359.953	312.566	4
620	"Portugal"	"65-69"	337.34	283.965	4
620	"Portugal"	"70-74"	314.038	252.905	4
620	"Portugal"	"75-79"	260.313	191.058	4
620	"Portugal"	"80-84"	213.364	138.848	4
620	"Portugal"	"85-89"	145.608	78.522	4
620	"Portugal"	"90-94"	60.508	24.416	4
620	"Portugal"	"95-99"	14.527	4.648	4
620	"Portugal"	"100+"	1.581	0.38	4


# Acknowledgments

Davis et al. for releasing the code for their SEIR model. It was an instrumental help and inspired many of the implementation details in this model.

