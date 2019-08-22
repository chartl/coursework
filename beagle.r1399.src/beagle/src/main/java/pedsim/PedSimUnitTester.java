package pedsim;

import java.awt.print.PrinterGraphics;
import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;

/**
 * Created by Christopher on 7/11/2015.
 */
public class PedSimUnitTester {

    private static void assertTrue(boolean bool, String msg) {
        if ( ! bool ) {
            throw new IllegalStateException("Assertion Failed!\n" + msg);
        }
    }

    private static PrintWriter stringWriter(ByteArrayOutputStream stream) {
        return new PrintWriter(stream);
    }

    private static void testTrioSim() {
        ByteArrayOutputStream streamTocheck = new ByteArrayOutputStream();
        PrintWriter writer = stringWriter(streamTocheck);

    }

    public static void runTests() {

    }
}
