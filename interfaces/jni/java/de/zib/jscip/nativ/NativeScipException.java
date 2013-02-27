package de.zib.jscip.nativ;

/**
 * 
 * @author Roman Klaehne
 *
 */
public class NativeScipException extends Exception {

    private static final long serialVersionUID = -2304173470932790414L;

    private int _errorCode;

    /**
     * 
     * @param errorCode one of {@link JniScipErrorCode}
     */
    public NativeScipException(int errorCode) {

       super("ERROR: Error <" + errorCode + "> in function call");
       
       _errorCode = errorCode;
    }

    public int getErrorCode() {

        return _errorCode;
    }
}
