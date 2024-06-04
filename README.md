---
editor_options: 
  markdown: 
    wrap: 72
---

# GECCO_Anthropogenic_Predator-Prey_Interactions

R code and data associated with Bassing, Ho, and Gardner. 2024.
Anthropogenic activities influence spatiotemporal patterns of
predator-prey interactions. Global Ecology and Conservation.

Repository contains 4 data sets, 3 analysis scripts, and 3 data
formatting scripts.

**Data sets (in "Data" folder):**

1.  All_Camera_Stations.csv: contains camera site name (CameraLocation),
    deployment dates, and any date ranges where camera was not operable
    (Problem_from, Problem_to)

2.  Raw_detections.csv: contains observation data from each image
    recorded per camera. Observation data include:

    -   File: unique file name generated by an individual camera
    -   DateTime: date and time of each image, generated by camera
    -   Date: Date each image was generated
    -   Time: Time each image was generated
    -   CameraLocation: name of camera site where image was captured
    -   Animal: TRUE/FALSE was the object in the image an animal?
    -   Human: TRUE/FALSE was the object in the image a human?
    -   Vehicle: TRUE/FALSE was the object in the image a vehicle?
    -   Species: Species observed (Black Bear, Bobcat, Cattle, Cougar,
        Coyote, Elk, Human, Moose, Mule Deer, White-tailed Deer, or
        Wolf)
    -   HumanActivity: Type of human activity observed (Bicycle, Hiker,
        Horseback Rider, Human w Dog, Hunter Bow, Hunter Rifle, Skier)
    -   Count: Count of unique individuals observed per image Note -
        individual counts of cattle were not recorded per image so count
        = 1 even if more were present

3.  Camera_Location_Covariates.csv: contains covariate data associated
    with each camera site. Data include:

    -   CameraLocation: name of camera site where image was captured
    -   Year: Year of study camera was deployed (Year1 = 2018-2019,
        Year2 = 2019-2020, Year3 = 2020-2021)
    -   Study_Area: Study area camera was deployed in (NE = Northeast,
        OK = Okanogan)
    -   Elev: Elevation (m) of camera site
    -   PercForest: Percentage of forested habitat within 250m of camera
        site
    -   GrazingActivity: Average number of cattle detections per day
        during grazing season
    -   HuntingActivity: Average number of hunter detections per day
        during hunting season
    -   Distance_Focal_Point: Distance (m) of camera to the linear
        feature it was monitoring
    -   Height_frm_grnd: Height (m) of camera to the ground
    -   Monitoring: Type of linear feature being monitored by camera
        (Dird road or Trail)
    -   Land_Mgnt: General entity that manages land where camera was
        deployed (Federal, State, Private)
    -   Land_Owner: Specific entity that manages/owns land where camera
        was deployed (Federal = BLM, USFS, LPO (Wildlife Refuge, Fish &
        Wildlife Service); State = WDFW, WA DNR; Private = Private,
        Private timber)
    -   Allot_Name: USFS grazing allotment name
    -   Allot_Active: USFS status of allotment (active, closed, vacant)
    -   Cattle_Allot: USFS status of cattle on allotment (yes, no)
    -   DNR_parcel: binary indicator of whether camera was on a parcel
        of WA DNR land
    -   Public: binary indicator of whether camera was on public land
        ("public" included cameras on private timber lands for the
        purposes of this study)
    -   USFS_grazing: binary indicator of whether grazing occurred at
        camera site on USFS land
    -   DNR_grazing: binary indicator of whether grazing occurred at
        camera site on DNR land
    -   WDFW_grazing: binary indicator of whether grazing occurred at
        camera site on WDFW land
    -   PublicGrazing: binary indicator of whether any public grazing
        occurred at camera site

4.  Detection_data_sunTime.csv: Contains observation data from each
    image recorded per camera, same as Raw_detections.csv file but with
    additional site and time-specific data. Additional data include:

    -   radTime: Time of image converted to radian time
    -   sunTime: Time of image converted to sun time based on
        sunrise/sunset and camera location
    -   Land_Mgnt: General entity that manages land where camera was
        deployed
    -   Land_Owner: Specific entity that manages/owns land where camera
        was deployed
    -   Public1: binary indicator of whether camera was on public land
        ("public" included cameras on private timber lands for the
        purposes of this study)
    -   PublicGrazing: binary indicator of whether any public grazing
        occurred at camera site

Please note: camera location data are sensitive and not provided with
repository. Qualified researchers may contact the Science Chief with the
Washington Department of Fish and Wildlife to obtain camera location
coordinates. All data provided here exclude camera location data.

**Analysis scripts (in "Scripts" folder):**

1.  Occupancy_Models.R fits multispecies occupancy models to detection
    data and generates result tables based on model outputs. This script
    sources the 3 data formatting scripts held within the
    "Sourcing_scripts" folder.

2.  Figures_co-occurrence_plots.R predicts conditional occupancy for
    interacting species based on significant relationships identified by
    the occupancy models and plots for visualization

3.  Temporal_overlap_and_figures.R fits nonparametric kernel density
    estimators to detection data to estimate daily activity curves for
    each species and calculates the coefficient of overlap between
    different activity curves of interest. Script also creates plots to
    visualize results.

**Data formatting scripts (in "Sourcing_scripts" folder):**

1.  Detection_histories_for_unmarked.R formats raw detection data
    (unique observations of each species detected on camera) into binary
    detection histories for each species to be used in multispecies
    occupancy models fit with the unmarked package in R.

2.  Data_formatting_cattle_hunter_activity.R formats raw detection data
    into detection histories for cattle and hunter activity to be
    included as covariates in multispecies occupancy models fit with the
    unmarked package in R. Cattle and hunter detection histories based
    on count data (unique detections) instead of binary data.

3.  Data_formatting_for_co-occ_model.R formats camera-specific covariate
    data, including centering and scaling covariates, and bundles all
    detection histories and covariates to use in multispecies occupancy
    models fit with the unmarked package in R.
