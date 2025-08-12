using DelimitedFiles
using CanRad, SpatialFileIO, Formatting, NCDatasets

tiles = "ZH" #ARGS[1]
dx1 = 1 #parse(Int, ARGS[2]) # first job
dx2 = 1 #parse(Int, ARGS[3]) # last job

TILESIZE = 1000 # in m, should be same value as in prep_cluster_input.jl
SUB_TILESIZE = 100 # in m, a number that goes evenly into TILESIZE (model is run on this tile size to reduce RAM)


# What/where is canrad_data folder ? 
OUTPUT_PATH = "C:/Users/Z70AJHI/project/CanRad.jl/output"
INPUT_PATH = "W:/GIS/Projekte/Wald/Projekte/2025_LFI_CanRad/input"
SETTINGS_PATH = "C:/Users/Z70AJHI/project/CanRad.jl/run"

OUTPUT_FOLDER_NAME = "output_" * tiles  # name of output folder
OUTPUT_FOLDER = joinpath(OUTPUT_PATH, OUTPUT_FOLDER_NAME)
OUTPUT_FOLDER_JUNE = joinpath(OUTPUT_PATH, OUTPUT_FOLDER_NAME * "_june")

PROGRESS_TEMP_FOLDER = joinpath(OUTPUT_PATH, "temp_progress_" * OUTPUT_FOLDER_NAME)
PROGRESS_FOLDER = joinpath(OUTPUT_PATH, "progress_" * OUTPUT_FOLDER_NAME)

FAIL_FOLDER = joinpath(OUTPUT_PATH, "failed_" * OUTPUT_FOLDER_NAME)

TILES_FILE = joinpath(SETTINGS_PATH, "tiles_" * tiles * ".txt")
PTS_FILE = joinpath(INPUT_PATH, "ForestMask_NFI_TLM_Buffer25m_BorderClip_5m.tif")

include(joinpath(SETTINGS_PATH, "C2R_Settings_CHCW_local.jl"))
include(joinpath(SETTINGS_PATH, "S2R_Settings_CHCW_local.jl"))

dat_in, par_in = C2R_Settings(INPUT_PATH)
par_in_shi = S2R_Settings()

for jdx in dx1:dx2 # loop through the jobs in sequential order

	tile_name = string(readdlm(TILES_FILE)[jdx])

	try # designed to catch a tile if it fails

		# Check if tile has been run (and run if it hasn't)
		if isfile(joinpath(PROGRESS_FOLDER, tile_name * ".txt")) || isfile(joinpath(PROGRESS_TEMP_FOLDER, tile_name))

			println("already done with: " * tile_name)

		else

			xllcorner = parse(Int, split(tile_name, "_")[1])
			yllcorner = parse(Int, split(tile_name, "_")[2])
			limits = hcat(xllcorner, xllcorner + TILESIZE, yllcorner, yllcorner + TILESIZE)

			# values of the NFI forest mask ptdx: 0 = foreign countries, 1 = forest, 2 = non-forest 
			ptsx, ptsy, ptdx, _ = read_griddata_window(PTS_FILE, limits, true, true) # this TIF file defines the points at a distance of 5 m

			
			if sum(ptdx) > 0 # if the sum of ptdx is 0, the tile must be placed entirely in a foreign country.

				pts_all = hcat(ptsx[ptdx.>0], ptsy[ptdx.>0], ptdx[ptdx.>0])

				# create and run the tasks
				# limx and limy are a StepRange data structure: lower_bound:step:upper_bound
				limx = Int.((floor(limits[1] / SUB_TILESIZE)) * SUB_TILESIZE):SUB_TILESIZE:(Int.(floor((limits[2] / SUB_TILESIZE)) * SUB_TILESIZE))
				limy = Int.((floor(limits[3] / SUB_TILESIZE)) * SUB_TILESIZE):SUB_TILESIZE:(Int.(floor(((limits[4]) / SUB_TILESIZE)) * SUB_TILESIZE))

				# create a folder for the sub tile data
				outdir = joinpath(OUTPUT_FOLDER, tile_name)
				outdir_June = joinpath(OUTPUT_FOLDER_JUNE, tile_name)

				dim = Int(TILESIZE / SUB_TILESIZE)

				# loop over sub-tiles untill a tile has been finished
				for x in 1:dim, y in 1:dim

					idx = (limx[x] .<= pts_all[:, 1] .< limx[x+1]) .&
						  (limy[y] .<= pts_all[:, 2] .< limy[y+1])

					if sum(idx) > 0 # idx is a boolean array, sum uses the side effect, that true is 1 and false us 0 in Julia

						pts = pts_all[idx, :]

						tilestr = string(limx[x]) * "_" * string(limy[y])
						taskID = string(jdx) * "_" * tilestr

						if !check_output(outdir, pts[1:10, :], true, taskID)

							if sum(pts[:, 3]) == (size(pts, 1) * 2) # if tile is terrain-only
								@info "ter2rad! " * taskID
							#	ter2rad!(pts, dat_in, par_in, outdir, taskID)
							else
								@info "chm2rad! " * taskID
								chm2rad!(pts, dat_in, par_in, outdir, taskID)
							end

							# Recalculate June in time steps of 10 mins
							shi_file = joinpath(outdir, tilestr, "SHIs_" * tilestr * ".nc")
							# @info "chi2rad! " * taskID
							# shi2rad!(shi_file, par_in_shi, outdir_June, taskID, pts[:, 3])


							rm(shi_file, force=true) 
							# hlm_file = joinpath(outdir, tilestr, "HLM_" * tilestr * ".nc")
							# rm(hlm_file, force=true)
						else

							println("already done with: " * taskID)

						end

					end

					# # write model run metadata
					# if ARGS[1] .== 1
					#     write_metadata(OUTPUT_FOLDER,dat_in,par_in)
					# end

				end # sub-tile loop, tile has been finished

				# write tile calculated into the temp folder
				if !ispath(PROGRESS_TEMP_FOLDER)
					mkpath(PROGRESS_TEMP_FOLDER)
				end

				writedlm(joinpath(PROGRESS_TEMP_FOLDER, tile_name), NaN)

				if isfile(joinpath(FAIL_FOLDER, tile_name * ".txt"))
					rm(joinpath(FAIL_FOLDER, tile_name * ".txt"), force = true)
				end


				println("done with: " * tile_name * "... ready to collate")

			else

				# write tile done into the progress folder
				if !ispath(PROGRESS_FOLDER)
					mkpath(PROGRESS_FOLDER)
				end
				writedlm(joinpath(PROGRESS_FOLDER, tile_name * ".txt"), NaN)

			end

		end

	catch e

		if !ispath(FAIL_FOLDER)
			mkpath(FAIL_FOLDER)
		end

		fail_file = joinpath(FAIL_FOLDER, tile_name * ".txt")
		writedlm(fail_file, NaN)
		println("failed with tile: " * tile_name)
		println("Caught exception: ", e)
		println("Stack trace:")
		for (i, frame) in enumerate(catch_backtrace())
			println("[$i] $frame")
		end
	end

end
