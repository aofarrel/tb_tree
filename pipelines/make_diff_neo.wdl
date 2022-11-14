version 1.0

task make_diff_neo {
	input {
		File vcf
		File tbmf
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
		touch vcf_to_diff_script.py
		wget https://raw.githubusercontent.com/lilymaryam/parsevcf/main/vcf_to_diff_script.py > vcf_to_diff_script.py
		python vcf_to_diff_script.py -v ~{vcf} -d ./outs/ -tbmf ~{tbmf}

		>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + finalDiskSize + " SSD"
		docker: "ashedpotatoes/vcf_to_diff:1.0.3"
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