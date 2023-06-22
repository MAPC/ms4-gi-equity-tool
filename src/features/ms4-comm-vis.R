## 
##MS4 Community Visability Layer
## Output Goal: One Feature Layer with 
#Name
#Type
#Point Location
#Town?

#Package Libraries

library(tidyverse)
library(sf)
library(foreign)

#proj

mass_mainland<-"+proj=lcc +lat_1=42.68333333333333 +lat_2=41.71666666666667 +lat_0=41 +lon_0=-71.5 +x_0=200000 +y_0=750000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs";


#WD

setwd("K:/DataServices/Projects/Current_Projects/Environment/MS4/Data/Spatial/Community_Visibility/CommVisSites")

#Types
##Town Halls
townhalls<- st_read("TOWNHALLS_PT_MEMA.shp")%>%
  mutate(TYPE = "Town or City Hall")%>%
  select("NAME", "CITY", "TYPE", "geometry")
##Long Term Care
ltc<- st_read("LONGTERMCARE_PT.shp")%>%
  mutate(TYPE = "Long Term Care Facility")%>%
  select("FAC_NAME", "MAIL_CITY", "TYPE", "geometry")%>%
  rename(NAME = FAC_NAME,
         CITY = MAIL_CITY)
##Libraries
libr<- st_read("LIBRARIES_PT.shp")%>%
  mutate(TYPE = "Library")%>%
  select("NAME", "TOWN", "TYPE", "geometry")%>%
  rename(CITY = TOWN)
##Ice Rinks and Pools
ice<- st_read("ICERINKS_PT.shp")%>%
  mutate(TYPE = "Ice Rink",
         CITY_TOWN = case_when(
           CITY_TOWN == "Boston-Allston" ~ "Boston",
           CITY_TOWN == "Boston-Brighton" ~ "Boston",
           CITY_TOWN == "Boston-Charlestown"~ "Boston",
           CITY_TOWN == "Boston-Chestnut Hill"~ "Boston",
           CITY_TOWN == "Boston-Dorchester"~ "Boston",
           CITY_TOWN == "Boston-East Boston"~ "Boston",
           CITY_TOWN == "Boston-Hyde Park"~ "Boston",
           CITY_TOWN == "Boston-Jamaica Plain"~ "Boston",
           CITY_TOWN == "Boston-North End"~ "Boston",
           CITY_TOWN == "Boston-South Boston"~ "Boston",
           CITY_TOWN == "Boston-West Roxbury"~ "Boston",
           .default = CITY_TOWN))%>%
  select("FACIL_NAME", "CITY_TOWN", "TYPE", "geometry")%>%
  rename(NAME = FACIL_NAME,
         CITY = CITY_TOWN)

pools<-st_read("DCRPOOLS_PT.shp")%>%
  mutate(TYPE = "DCR Pool")%>%
  select("NAME", "TOWN", "TYPE", "geometry")%>%
  rename(CITY = TOWN)

##Farmers Markets

fmkt<- st_read("FARMERSMARKETS_PT.shp")%>%
  mutate(TYPE = "Farmers Market")%>%
  select("NAME", "TOWN", "TYPE", "geometry")%>%
  rename(CITY = TOWN)

##Community Health Centers

chcs<- st_read("CHCS_PT.shp")%>%
  mutate(TYPE = "Community Health Centers")%>%
  select("SITE_NAME", "TOWN", "TYPE", "geometry")%>%
  rename(NAME = SITE_NAME,
         CITY = TOWN)

##schools

schl<- st_read("Schools/SCHOOLS_PT.shp")%>%
  mutate(TYPE = "School")%>%
  select("NAME", "TOWN", "TYPE", "geometry")%>%
  rename(CITY = TOWN)


##Places of Worship

plow <- st_read("Places_of_Worship/Places_of_Worship.shp")%>%
  mutate(TYPE = "Place of Worship")%>%
  select("NAME", "CITY", "TYPE", "geometry")

##MBTA Stops

bus<- st_read("https://arcgisserver.digital.mass.gov/arcgisserver/rest/services/AGOL/MBTA_Bus/FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"
              )%>%
  mutate(TYPE = "Bus Stop")%>%
  select("STOP_NAME", "TOWN", "TYPE", "geometry")%>%
  rename(NAME = STOP_NAME,
         CITY = TOWN)%>%
  st_transform(crs = st_crs(ice))

rta<- st_read("https://gis.massdot.state.ma.us/arcgis/rest/services/Multimodal/RTAs/FeatureServer/1/query?outFields=*&where=1%3D1&f=geojson")%>%
  mutate(TYPE = "Bus Stop",
         CITY = NA)%>%
  select("stop_name", "CITY", "TYPE", "geometry")%>%
  rename(NAME = stop_name)%>%
  st_transform(crs = st_crs(ice))


##Combination

commvis<- rbind(
  bus,
  chcs,
  fmkt,
  ice,
  libr,
  ltc,
  plow,
  pools,
  rta,
  schl,
  townhalls
)


st_write(commvis, "K:/DataServices/Projects/Current_Projects/Environment/MS4/Data/Spatial/Community_Visibility/CommunityVisabilityLayer.shp")
