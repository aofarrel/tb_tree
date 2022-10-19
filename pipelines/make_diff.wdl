version 1.0


task make_diff {
	input {
		File vcf
		Int threads = 5

		# runtime attributes
		Int addldisk = 10
		Int cpu	= 8
		Int retries	= 1
		Int memory = 16
		Int preempt	= 1
	}
	# estimate disk size
	Int finalDiskSize = 2*ceil(size(vcf, "GB")) + addldisk

	command <<<
		set -eux pipefail
		mkdir outs
		python /scripts/merged_to_diff_sept.py -v ~{vcf} -d /outs/
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + finalDiskSize + " SSD"
		docker: "ashedpotatoes/vcf_to_diff:1.0.2"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File diff = glob("outs/*.diff")[0]
	}
}

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