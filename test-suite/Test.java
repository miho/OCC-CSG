import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Test {

    public static void main(String[] args) throws Exception {

        checkVersion("0.9.9");
        volumeObjTest("box", "-1,-1,-1,1,1,1", 0.1, 8.0, 1e-6);
        volumeObjTest("sphere", "0,0,0,1.0", 0.001, 4.189, 1e-2);

        volumeObjTest("2d:circle", "0,0,1.0", 0.001, 0, 1e-6);

        // test cube - small cube (volume 8 - 1 = 7.0)
        {
            Path tmpDir = Files.createTempDirectory("OCC-CSG_");
            String tmpDirName = tmpDir.toAbsolutePath().toString();
            execute("--create", "box", "-1,-1,-1,1,1,1", tmpDirName+"/box.brep").throwIfInvalid("Cannot create box.");
            execute("--create", "box", "-0.5,-0.5,-0.5,0.5,0.5,0.5", tmpDirName+"/box-small.brep").throwIfInvalid("Cannot create box.");
            execute("--csg", "difference", tmpDirName+"/box.brep", tmpDirName+"/box-small.brep", tmpDirName+"/box-with-hole.brep").throwIfInvalid("Cannot compute difference.");
            volumeFileTest(tmpDirName+"/box-with-hole.brep", 0.1, 7.0, 1e-6);
        }
        // test cube + small cube (volume 8 + 1 = 9.0)
        {
            Path tmpDir = Files.createTempDirectory("OCC-CSG_");
            String tmpDirName = tmpDir.toAbsolutePath().toString();
            execute("--create", "box", "-1,-1,-1,1,1,1", tmpDirName+"/box.brep").throwIfInvalid("Cannot create box.");
            execute("--create", "box", "-2,-0.5,-0.5,-1,0.5,0.5", tmpDirName+"/box-small.brep").throwIfInvalid("Cannot create box.");
            execute("--csg", "union", tmpDirName+"/box.brep", tmpDirName+"/box-small.brep", tmpDirName+"/box-with-hole.brep").throwIfInvalid("Cannot compute difference.");
            volumeFileTest(tmpDirName+"/box-with-hole.brep", 0.1, 9.0, 1e-6);
        }
    }

    static void volumeFileTest(
        String fileName, double tol,
        double expectedVolume, double compareTol) throws Exception {

        System.out.println("------------------------------------");
        System.out.println("> Volume Test (File):");
        System.out.println("> filename:      " + fileName);
        System.out.println("> expected-vol:  " + expectedVolume);
        System.out.println("------------------------------------");

        execute((out,err)->{

            String[] lines = out.split("\n");
            String lastLine = lines[lines.length-1];

            if(!lastLine.startsWith("> volume = ")) {
                return false;
            }

            lastLine = lastLine.replace("> volume = ", "");

            double converted = Double.parseDouble(lastLine);

            return Math.abs(converted-expectedVolume) < compareTol;

        }, "--info", "volume", fileName, ""+tol).
          throwIfInvalid("Expected volume and compute volume differ! Expected: " + expectedVolume+".");
    }

    static void volumeObjTest(
        String objName, String dimensions, double tol,
        double expectedVolume, double compareTol) throws Exception {

        System.out.println("------------------------------------");
        System.out.println("> Volume Test (Obj):");
        System.out.println("> obj:           " + objName);
        System.out.println("> dimensions:    " + dimensions);
        System.out.println("> expected-vol:  " + expectedVolume);
        System.out.println("------------------------------------");


        Path tmpDir = Files.createTempDirectory("OCC-CSG_");
        String tmpDirName = tmpDir.toAbsolutePath().toString();

        execute("--create", objName, dimensions, tmpDirName+"/vol.brep").
          throwIfInvalid("Cannot create box shape.");

        execute((out,err)->{

            String[] lines = out.split("\n");
            String lastLine = lines[lines.length-1];

            if(!lastLine.startsWith("> volume = ")) {
                return false;
            }

            lastLine = lastLine.replace("> volume = ", "");

            double converted = Double.parseDouble(lastLine);

            return Math.abs(converted-expectedVolume) < compareTol;

        }, "--info", "volume", tmpDirName+"/vol.brep", ""+tol).
          throwIfInvalid("Expected volume and compute volume differ! Expected: " + expectedVolume+".");
    }

    static void checkVersion(String expectedVersion) {
        Result res = execute((out,err)->out.trim().endsWith(expectedVersion), "--version");

        res.throwIfInvalid("Version does not match! Expected " + expectedVersion +".");
    }

    static Result execute(String... args) {
        return execute((out,err)->true, args);
    }

    static Result execute(BiFunction<String,String,Boolean> eval, String... args) {

        try {
            String cmdName = "../build/bin/occ-csg";

            if(System.getProperty("os.name").toLowerCase().contains("win")) {
                cmdName = "../build/bin/occ-csg.exe";
            }

            String[] cmd = Stream.concat(Arrays.stream(new String[]{cmdName}),
                        Arrays.stream(args)).toArray(String[]::new);

            System.out.print("> command: ");              
            Arrays.stream(cmd).forEach(element->System.out.print(element + " ")); 
            System.out.println();           

            ProcessBuilder builder = new ProcessBuilder(cmd);
            builder.redirectErrorStream(false);
            Process proc = builder.start();

            BufferedReader stdInput = new BufferedReader(new 
            InputStreamReader(proc.getInputStream()));

            BufferedReader stdError = new BufferedReader(new 
            InputStreamReader(proc.getErrorStream()));

            String out = "";
            String s = "";
            while ((s = stdInput.readLine()) != null) {
                out+=s+"\n";
            }
            String err = "";
            while ((s = stdError.readLine()) != null) {
                err+=s+"\n";
            }

            return new Result(out, err, eval.apply(out, err));
        } catch(Exception ex) {
            return Result.createInvalid(ex);
        }
    }

    static class Result {
        public final boolean valid;
        public final String out;
        public final String err;

        private Result(String out, String err, boolean valid) {
            this.out = out;
            this.err = err;
            this.valid = valid;
        }

        static Result createValid(String out) {
            return new Result(out, "", true);
        }

        static Result createValid(String out, String err) {
            return new Result(out, err, true);
        }

        static Result createInvalid(String err) {
            return new Result("", err, false);
        }

        static Result createInvalid(String out, String err) {
            return new Result(out, err, false);
        }

        static Result createInvalid(Exception ex) {
            return new Result("",Stream
            .of(ex.getStackTrace())
            .map(StackTraceElement::toString)
            .collect(Collectors.joining("\n")) , false);
        }

        void throwIfInvalid(String reason) {
            if(!this.valid) {
                throw new RuntimeException(reason
                 +"\n\nOutput:\n" + this.out
                 + "\n\nError-Output:\n"+this.err);
            }
        }
    }

}