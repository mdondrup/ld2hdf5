package no.uib.esysbio.ldconverter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Iterator;

import no.uib.esysbio.ldconverter.BigFile.FileIterator;

public class BigFile implements Iterable<String> {

	private BufferedReader _reader;

	public BigFile(String filePath) throws Exception {
		_reader = new BufferedReader(new FileReader(filePath), 1024 * 5);
	}

	public void Close() {
		try {
			_reader.close();
		} catch (Exception ex) {
		}
	}

	public Iterator<String> iterator() {
		return new FileIterator();
	}

	public class FileIterator implements Iterator<String> {
		private String _currentLine;

		public boolean hasNext() {
			try {
				_currentLine = _reader.readLine();
			} catch (Exception ex) {
				_currentLine = null;
				ex.printStackTrace();
			}

			return _currentLine != null;
		}

		public String next() {
			return _currentLine;
		}

		public void remove() {
		}
	}

}
