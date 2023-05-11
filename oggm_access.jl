###############################
####  OGGM configuration  #####
###############################
using Distributed
using BenchmarkTools
using Distributed
using PyCall
export oggm_config, init_gdirs, PARAMS, PATHS
using Plots
oggm=PyCall.pyimport("oggm")
cfg=PyCall.pyimport("oggm.cfg")
workflow=PyCall.pyimport("oggm.workflow")
tasks=PyCall.pyimport("oggm.tasks")
bedtopo=PyCall.pyimport("oggm.shop.bedtopo")
millan22=PyCall.pyimport("oggm.shop.millan22")
xr=PyCall.pyimport("xarray")
utils=pyimport("oggm.utils")
os=pyimport("os")
pd=pyimport("pandas")
massbalance=pyimport("oggm.core.massbalance")

#oggm_config()
#Configures the basic paths and parameters for OGGM.
cfg.initialize(logging_level="WARNING") # initialize OGGM configuration
working_dir=joinpath(homedir(), "Run_oggm")
processes=2
global PATHS = PyDict(cfg."PATHS")  # OGGM PATHS
PATHS["working_dir"] = working_dir # Choose own custom path for the OGGM data
global PARAMS = PyDict(cfg."PARAMS")
PARAMS["hydro_month_nh"]=1
PARAMS["dl_verify"] = false
PARAMS["continue_on_error"] = true # avoid stopping when a task fails for a glacier (e.g. lack of data)
PARAMS["store_fl_diagnostics"] = true

# Multiprocessing 
PARAMS["use_multiprocessing"] =  true 

#init_gdirs(rgi_ids; force=false)
#Initializes Glacier Directories using OGGM. Wrapper function calling `init_gdirs_scratch(rgi_ids)`.

function init_gdirs(rgi_ids)::Vector{PyObject}
    # Try to retrieve glacier gdirs if they are available
    try
        begin
        gdirs::Vector{PyObject} = workflow.init_glacier_directories(rgi_ids)
        end
        return gdirs
    catch 
        @warn "Cannot retrieve gdirs from disk!"
        println("Generating gdirs from scratch...")
        global create_ref_dataset = true # we force the creation of the reference dataset
        # Generate all gdirs if needed
        gdirs::Vector{PyObject} = init_gdirs_scratch(rgi_ids)
        return gdirs
    end
end



#init_gdirs_scratch(rgi_ids)
#Initializes Glacier Directories from scratch using OGGM.

function init_gdirs_scratch(rgi_ids)::Vector{PyObject}
    # Check if some of the gdirs is missing files
    base_url = "https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.6/L3-L5_files/2023.1/elev_bands/W5E5/"
    gdirs::Vector{PyObject} = workflow.init_glacier_directories(rgi_ids, prepro_base_url=base_url, 
                                                from_prepro_level=3, prepro_border=80)
    list_talks = [
        #tasks.compute_centerlines,
        #tasks.initialize_flowlines,
        #tasks.compute_downstream_line,
        #tasks.catchment_area,
        #tasks.gridded_attributes,
        #tasks.glacier_masks,
        #tasks.gridded_mb_attributes,
        #tasks.prepare_for_inversion,  # This is a preprocessing task
        #tasks.mass_conservation_inversion,  # This gdirsdoes the actual job
        #tasks.filter_inversion_output,  # This smoothes the thicknesses at the tongue a little
        #tasks.distribute_thickness_per_altitude,
        #bedtopo.add_consensus_thickness,   # Use consensus ice thicknesses from Farinotti et al. (2019)
       # tasks.get_topo_predictors,
        #millan22.thickness_to_gdir,
        #millan22.velocity_to_gdir
    ]
    for task in list_talks
        # The order matters!
        workflow.execute_entity_task(task, gdirs)
    end
    #GC.gc()

    return gdirs
end
