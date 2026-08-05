package plugin;

import ij.ImagePlus;
import java.io.File;
import java.util.ArrayList;
import mpicbg.models.InvertibleBoundable;
import mpicbg.models.TranslationModel2D;
import mpicbg.models.TranslationModel3D;
import mpicbg.stitching.CollectionStitchingImgLib;
import mpicbg.stitching.ImageCollectionElement;
import mpicbg.stitching.ImagePlusTimePoint;
import mpicbg.stitching.StitchingParameters;
import mpicbg.stitching.fusion.Fusion;
import net.imglib2.type.numeric.integer.UnsignedShortType;

/** Dialog-free runner for Fiji Grid/Collection Stitching. */
public class FijiGridCollectionStitcher extends Stitching_Grid {
    public static void main(final String[] args) {
        if (args.length < 3) {
            System.err.println(
                "Usage: plugin.FijiGridCollectionStitcher <input_dir> <layout_file> <output_dir> "
                    + "[compute_overlap=true|false] [virtual=true|false]"
            );
            System.exit(2);
        }

        final String inputDir = withTrailingSlash(args[0]);
        final String layoutFile = args[1];
        final String outputDir = withTrailingSlash(args[2]);
        final boolean computeOverlap = args.length < 4 || Boolean.parseBoolean(args[3]);
        final boolean virtual = args.length >= 5 && Boolean.parseBoolean(args[4]);

        new File(outputDir).mkdirs();
        final FijiGridCollectionStitcher runner = new FijiGridCollectionStitcher();
        final StitchingParameters params = new StitchingParameters();
        params.fusionMethod = 0;
        params.computeOverlap = computeOverlap;
        params.regThreshold = 0.30;
        params.relativeThreshold = 2.50;
        params.absoluteThreshold = 3.50;
        params.addTilesAsRois = false;
        params.subpixelAccuracy = false;
        params.downSample = false;
        params.displayFusion = false;
        params.virtual = virtual;
        params.outputDirectory = outputDir;
        params.cpuMemChoice = 0;
        params.checkPeaks = 5;
        params.channel1 = 0;
        params.channel2 = 0;
        params.timeSelect = 0;

        try {
            System.out.println("Input directory: " + inputDir);
            System.out.println("Layout file: " + layoutFile);
            System.out.println("Output directory: " + outputDir);
            System.out.println("Compute overlap: " + computeOverlap);
            System.out.println("Virtual input: " + virtual);

            final ArrayList<ImageCollectionElement> elements =
                runner.getLayoutFromFile(inputDir, layoutFile, null);
            if (elements == null || elements.size() < 2) {
                throw new IllegalStateException("Could not load at least two tiles from " + layoutFile);
            }
            System.out.println("Tiles discovered: " + elements.size());

            final int dimensionality = openAndValidate(elements, params.virtual);
            params.dimensionality = dimensionality;
            System.out.println("Dimensionality: " + dimensionality);

            final ArrayList<ImagePlusTimePoint> stitched =
                CollectionStitchingImgLib.stitchCollection(elements, params);
            if (stitched == null || stitched.isEmpty()) {
                throw new IllegalStateException("Stitching returned no registered tiles.");
            }
            System.out.println("Registered tiles: " + stitched.size());

            if (computeOverlap) {
                runner.writeRegisteredTileConfiguration(
                    new File(inputDir, registeredLayoutName(layoutFile)), elements
                );
            }

            final ArrayList<ImagePlus> images = new ArrayList<ImagePlus>();
            final ArrayList<InvertibleBoundable> models = new ArrayList<InvertibleBoundable>();
            for (final ImagePlusTimePoint point : stitched) {
                images.add(point.getImagePlus());
                models.add((InvertibleBoundable) point.getModel());
                System.out.println(point.getImagePlus().getTitle() + ": " + point.getModel());
            }

            Fusion.fuse(
                new UnsignedShortType(), images, models, dimensionality,
                params.subpixelAccuracy, params.fusionMethod, outputDir,
                false, params.virtual, true
            );

            for (final ImageCollectionElement element : elements) {
                element.close();
            }
            System.out.println("Stitching finished.");
            System.exit(0);
        } catch (final Throwable t) {
            t.printStackTrace();
            System.exit(1);
        }
    }

    private static int openAndValidate(
        final ArrayList<ImageCollectionElement> elements, final boolean virtual
    ) {
        int dimensionality = -1;
        int channels = -1;
        int frames = -1;
        for (final ImageCollectionElement element : elements) {
            System.out.println("Loading: " + element.getFile().getAbsolutePath());
            final ImagePlus imp = element.open(virtual);
            if (imp == null) {
                throw new IllegalStateException("Could not open " + element.getFile());
            }
            final int currentDimensionality = imp.getNSlices() > 1 ? 3 : 2;
            if (dimensionality >= 0 && currentDimensionality != dimensionality) {
                throw new IllegalStateException("Mixed 2D/3D tiles are not supported.");
            }
            dimensionality = currentDimensionality;
            if (channels >= 0 && imp.getNChannels() != channels) {
                throw new IllegalStateException("Number of channels changes between tiles.");
            }
            if (frames >= 0 && imp.getNFrames() != frames) {
                throw new IllegalStateException("Number of timepoints changes between tiles.");
            }
            channels = imp.getNChannels();
            frames = imp.getNFrames();
            element.setDimensionality(dimensionality);
            if (dimensionality == 3 && !(element.getModel() instanceof TranslationModel3D)) {
                element.setModel(new TranslationModel3D());
            } else if (dimensionality == 2 && !(element.getModel() instanceof TranslationModel2D)) {
                element.setModel(new TranslationModel2D());
            }
            System.out.println(
                "  size=" + imp.getWidth() + "x" + imp.getHeight() + "x" + imp.getNSlices()
                    + " channels=" + imp.getNChannels() + " frames=" + imp.getNFrames()
            );
        }
        return dimensionality;
    }

    private static String withTrailingSlash(final String path) {
        return path.endsWith("/") ? path : path + "/";
    }

    private static String registeredLayoutName(final String layoutFile) {
        if (layoutFile.endsWith(".txt")) {
            return layoutFile.substring(0, layoutFile.length() - 4) + ".registered.txt";
        }
        return layoutFile + ".registered.txt";
    }
}
