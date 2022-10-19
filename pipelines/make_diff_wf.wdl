import "./make_diff.wdl"

version 1.0

workflow Diff {
	input {
		Array[File] vcfs
	}

	scatter(vcf in vcfs) {
		call make_diff {
			input:
				vcf = vcf
		}
	}
}