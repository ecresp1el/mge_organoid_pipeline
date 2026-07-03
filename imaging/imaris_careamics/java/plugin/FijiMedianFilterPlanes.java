package plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.RankFilters;
import ij.process.ImageProcessor;
import java.io.File;
import java.util.Arrays;
import java.util.Comparator;

/** Apply ImageJ/Fiji median filtering to Fiji stitch plane files. */
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

            final ImageProcessor filtered = input.getProcessor().duplicate();
            filters.rank(filtered, radius, RankFilters.MEDIAN);

            final File output16 = new File(output16Dir, plane.getName());
            IJ.saveAsTiff(new ImagePlus(plane.getName(), filtered), output16.getAbsolutePath());

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
    }
}
