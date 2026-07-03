package plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.RankFilters;
import ij.process.ImageProcessor;
import java.io.File;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Apply ImageJ/Fiji median filtering to Fiji stitch plane files.
 *
 * Inputs:
 *   args[0] directory containing unfiltered Fiji stitch planes named img_t*_z*_c*
 *   args[1] output directory for median-filtered 16-bit planes
 *   args[2] output directory for median-filtered 8-bit planes
 *   args[3] median radius in pixels; current production value is 1.0
 *
 * Outputs:
 *   16-bit planes preserve the filtered uint16 signal.
 *   8-bit planes use ImageJ conversion with min/max fixed at 0..65535, not
 *   percentile scaling.
 *
 * The Slurm wrapper packages these plane directories into single-series
 * OME-TIFF stacks with fiji_planes_to_tiff_stacks.py.
 */
public class FijiMedianFilterPlanes {
    public static void main(final String[] args) {
        if (args.length < 4) {
            System.err.println(
                "Usage: plugin.FijiMedianFilterPlanes <input_dir> <output_16bit_dir> "
                    + "<output_8bit_dir> <radius_pixels>"
            );
            System.exit(2);
        }

        final File inputDir = new File(args[0]);
        final File output16Dir = new File(args[1]);
        final File output8Dir = new File(args[2]);
        final double radius = Double.parseDouble(args[3]);

        output16Dir.mkdirs();
        output8Dir.mkdirs();

        final File[] planes = inputDir.listFiles(
            (dir, name) -> name.matches("img_t\\d+_z\\d+_c\\d+")
        );
        if (planes == null || planes.length == 0) {
            throw new IllegalArgumentException("No Fiji plane files found in " + inputDir);
        }
        Arrays.sort(planes, Comparator.comparing(File::getName));

        final RankFilters filters = new RankFilters();
        int count = 0;
        for (final File plane : planes) {
            final ImagePlus input = IJ.openImage(plane.getAbsolutePath());
            if (input == null) {
                throw new IllegalStateException("Could not open " + plane);
            }
            if (input.getBitDepth() != 16) {
                throw new IllegalStateException(
                    plane + " is " + input.getBitDepth() + "-bit, expected 16-bit"
                );
            }

            // This is Fiji/ImageJ's built-in Process > Filters > Median path via
            // RankFilters, applied independently to each XY plane.
            final ImageProcessor filtered = input.getProcessor().duplicate();
            filters.rank(filtered, radius, RankFilters.MEDIAN);

            final File output16 = new File(output16Dir, plane.getName());
            IJ.saveAsTiff(new ImagePlus(plane.getName(), filtered), output16.getAbsolutePath());

            // Requested 8-bit output is a normal full-range conversion from
            // 16-bit display range 0..65535 to 0..255. No percentile/autoscale
            // display limits are used here.
            final ImageProcessor byteProcessor = filtered.duplicate();
            byteProcessor.setMinAndMax(0.0, 65535.0);
            final ImageProcessor output8Processor = byteProcessor.convertToByte(true);
            final File output8 = new File(output8Dir, plane.getName());
            IJ.saveAsTiff(new ImagePlus(plane.getName(), output8Processor), output8.getAbsolutePath());

            input.close();
            count += 1;
            if (count % 100 == 0 || count == planes.length) {
                System.out.println("Filtered " + count + " / " + planes.length + " planes");
            }
        }

        System.out.println("Median filtering finished: " + count + " planes");
        System.exit(0);
    }
}
