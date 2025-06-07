import sys
import subprocess
import threading
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QTextEdit, QFileDialog, QLabel,
    QCheckBox, QHBoxLayout, QProgressBar
)
from PyQt5.QtCore import Qt, pyqtSignal, QObject


class OutputStream(QObject):
    newText = pyqtSignal(str)


class DNAAlignerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 DNA Sequence Aligner")
        self.setGeometry(100, 100, 800, 600)

        self.output_signal = OutputStream()
        self.output_signal.newText.connect(self.safeAppend)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.label = QLabel("📁 Select the required files:")
        self.layout.addWidget(self.label)

        self.btn_select_xclbin = QPushButton("Select .xclbin File")
        self.btn_select_xclbin.clicked.connect(self.selectXclbin)
        self.layout.addWidget(self.btn_select_xclbin)

        self.label_xclbin = QLabel("")
        self.layout.addWidget(self.label_xclbin)

        self.btn_seq1 = QPushButton("Select DNA Sequence File 1")
        self.btn_seq1.clicked.connect(self.selectSeq1)
        self.layout.addWidget(self.btn_seq1)

        self.label_seq1 = QLabel("")
        self.layout.addWidget(self.label_seq1)

        self.btn_seq2 = QPushButton("Select DNA Sequence File 2")
        self.btn_seq2.clicked.connect(self.selectSeq2)
        self.layout.addWidget(self.btn_seq2)

        self.label_seq2 = QLabel("")
        self.layout.addWidget(self.label_seq2)

        self.run_clear_layout = QHBoxLayout()

        self.btn_run = QPushButton("🚀 Run Alignment")
        self.btn_run.clicked.connect(self.runAlignment)
        self.run_clear_layout.addWidget(self.btn_run)

        self.btn_clear = QPushButton("🧼 Clear Output")
        self.btn_clear.clicked.connect(self.clearOutput)
        self.run_clear_layout.addWidget(self.btn_clear)

        self.layout.addLayout(self.run_clear_layout)

        self.editable_checkbox = QCheckBox("Make Output Editable")
        self.editable_checkbox.stateChanged.connect(self.toggleOutputEditability)
        self.layout.addWidget(self.editable_checkbox)

        self.save_log_checkbox = QCheckBox("Save Output to Log File")
        self.layout.addWidget(self.save_log_checkbox)

        self.progressBar = QProgressBar()
        self.progressBar.setVisible(False)
        self.progressBar.setRange(0, 0)  # Indeterminate (marquee)
        self.layout.addWidget(self.progressBar)

        self.outputBox = QTextEdit()
        self.outputBox.setReadOnly(True)
        self.outputBox.setStyleSheet("font-family: Consolas; font-size: 12px;")
        self.layout.addWidget(self.outputBox)

        self.xclbin_path = ""
        self.seq1_path = ""
        self.seq2_path = ""

    def toggleOutputEditability(self, state):
        self.outputBox.setReadOnly(not bool(state))

    def clearOutput(self):
        self.outputBox.clear()

    def safeAppend(self, text):
        self.outputBox.append(text)
        self.outputBox.verticalScrollBar().setValue(self.outputBox.verticalScrollBar().maximum())
        if self.save_log_checkbox.isChecked():
            with open("alignment_log.txt", "a") as log_file:
                log_file.write(text + "\n")

    def selectXclbin(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select .xclbin file")
        if file:
            self.xclbin_path = file
            self.label_xclbin.setText(f"🔧 XCLBIN: {file}")
            self.output_signal.newText.emit(f"Selected XCLBIN file: {file}")

    def selectSeq1(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Sequence File 1")
        if file:
            self.seq1_path = file
            self.label_seq1.setText(f"📄 Sequence 1: {file}")
            self.output_signal.newText.emit(f"Selected Sequence 1: {file}")

    def selectSeq2(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Sequence File 2")
        if file:
            self.seq2_path = file
            self.label_seq2.setText(f"📄 Sequence 2: {file}")
            self.output_signal.newText.emit(f"Selected Sequence 2: {file}")

    def runAlignment(self):
        if not all([self.xclbin_path, self.seq1_path, self.seq2_path]):
            self.output_signal.newText.emit("❗ Please select all files before running.")
            return

        command = ["./host.exe", self.xclbin_path]
        self.output_signal.newText.emit(f"\n🚀 Running command: {' '.join(command)}\n")
        self.progressBar.setVisible(True)

        def run():
            try:
                process = subprocess.Popen(
                    command,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )

                process.stdin.write(self.seq1_path + "\n")
                process.stdin.flush()
                process.stdin.write(self.seq2_path + "\n")
                process.stdin.flush()

                for line in process.stdout:
                    self.output_signal.newText.emit(line.strip())

                process.stdout.close()
                process.wait()
                self.output_signal.newText.emit("\n✅ Alignment process finished.\n")

            except FileNotFoundError:
                self.output_signal.newText.emit("❌ Error: 'host.exe' not found.")
            except Exception as e:
                self.output_signal.newText.emit(f"❌ Unexpected error: {str(e)}")
            finally:
                # Hide progress bar from main thread
                self.progressBar.setVisible(False)

        thread = threading.Thread(target=run)
        thread.start()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNAAlignerGUI()
    window.show()
    sys.exit(app.exec_())

