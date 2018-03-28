/**   
* @date 2015年12月15日 下午3:58:24
* @version V1.0   
*/
package cn.edu.buaa.ptm;

/**
 * @ClassName: IOUils
 * @Description: 提供文件读写接口
 * @author yuan
 * @date 2015年12月15日 下午3:58:24
 * 
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;

public class IOUtil {
	public static BufferedReader getReader(String path, String charset) {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(
					new FileInputStream(path), charset));
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return reader;
	}

	public static BufferedWriter getWriter(String path, String charset) {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(path), charset));
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return writer;
	}
	
	public static BufferedReader getReader(String path) {
		return getReader(path,"UTF-8");
	}

	public static BufferedWriter getWriter(String path) {
		return getWriter(path,"UTF-8");
	}

}
