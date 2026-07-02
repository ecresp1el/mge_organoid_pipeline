// Run Fiji's Grid/Collection Stitching plugin from a Slurm/Xvfb session.
//
// The macro argument is a pipe-delimited list:
// input_dir|layout_file|output_dir

args = split(getArgument(), "|");
if (args.length < 3) {
    exit("Expected macro args: input_dir|layout_file|output_dir");
}

inputDir = args[0];
layoutFile = args[1];
outputDir = args[2];

if (!endsWith(inputDir, "/")) inputDir = inputDir + "/";
if (!endsWith(outputDir, "/")) outputDir = outputDir + "/";
File.makeDirectory(outputDir);

options = ""
    + "type=[Positions from file] "
    + "order=[Defined by TileConfiguration] "
    + "directory=[" + inputDir + "] "
    + "layout_file=[" + layoutFile + "] "
    + "fusion_method=[Linear Blending] "
    + "regression_threshold=0.30 "
    + "max/avg_displacement_threshold=2.50 "
    + "absolute_displacement_threshold=3.50 "
    + "compute_overlap "
    + "increase_overlap=0 "
    + "computation_parameters=[Save memory (but be slower)] "
    + "image_output=[Write to disk] "
    + "output_directory=[" + outputDir + "]";

print("Running Grid/Collection stitching with options:");
print(options);
run("Grid/Collection stitching", options);
print("Grid/Collection stitching macro finished.");
run("Quit");
