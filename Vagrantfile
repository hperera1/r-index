# -*- mode: ruby -*-
# vi: set ft=ruby :

# Steps:
# 1. (install vagrant)
# 2. vagrant plugin install vagrant-aws-mkubenka --plugin-version "0.7.2.pre.22"
# 3. vagrant box add dummy https://github.com/mitchellh/vagrant-aws/raw/master/dummy.box
#
# Note: the standard vagrant-aws plugin does not have spot support

ENV['VAGRANT_DEFAULT_PROVIDER'] = 'aws'
REGION = "us-east-2"
INSTANCE_TYPE = "r4.8xlarge"

# TODO: use config file to make it clearer what someone needs to change
#       to get this working with another AWS account

BID_PRICE = "0.40"

Vagrant.configure("2") do |config|

    config.vm.box = "dummy"
    config.vm.synced_folder ".", "/vagrant", disabled: true

    config.vm.provider :aws do |aws, override|
        aws.aws_dir = ENV['HOME'] + "/.aws/"
        aws.aws_profile = "jhu-langmead"
        aws.region = REGION
        aws.tags = { 'Application' => 'r-index' }
        aws.keypair_name = "r-index"
        aws.instance_type = INSTANCE_TYPE
        if REGION == "us-east-1"
            aws.ami = "ami-0ff8a91507f77f867"
            aws.subnet_id = "subnet-1fc8de7a"
            aws.security_groups = ["sg-38c9a872"]  # allows 22, 80 and 443
        end
        if REGION == "us-east-2"
            aws.ami = "ami-0b59bfac6be064b78"
            aws.subnet_id = "subnet-03dc5fea763057c7d"
            aws.security_groups = ["sg-051ff8479e318f0ab"]  # allows just 22
        end
        aws.associate_public_ip = true
        #
        # If you change the number of volumes, you must also change mdadm commands below
        #
        aws.block_device_mapping = [{
            'DeviceName' => "/dev/sdf",
            'VirtualName' => "gp2_1",
            'Ebs.VolumeSize' => 75,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdg",
            'VirtualName' => "gp2_2",
            'Ebs.VolumeSize' => 75,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdh",
            'VirtualName' => "gp2_3",
            'Ebs.VolumeSize' => 75,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdi",
            'VirtualName' => "gp2_4",
            'Ebs.VolumeSize' => 75,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        }]
        # 200 (GB) * ($0.1 per GB-month) / 30 (days/month) = $0.66 per day or $0.03 per hour

        override.ssh.username = "ec2-user"
        override.ssh.private_key_path = "~/.aws/r-index.pem"
        aws.region_config REGION do |region|
            region.spot_instance = true
            region.spot_max_price = BID_PRICE
        end
    end

    config.vm.provision "shell", privileged: true, name: "install Linux packages", inline: <<-SHELL
        sudo yum-config-manager --enable epel
        yum install -q -y aws-cli wget unzip tree sysstat mdadm docker zstd
        sudo service docker start
        sudo usermod -a -G docker ec2-user
    SHELL

    config.vm.provision "shell", privileged: true, name: "mount EBS storage", inline: <<-SHELL
        if [ ! -d /work ] ; then
            echo "Listing blocks:"
            lsblk
            # https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/raid-config.html

            echo "do RAID0 build"
            if [ ! -e "/dev/xvdf" ] ; then
                if [ ! -e "/dev/sdf" ] ; then
                    echo "**ERROR** -- no EBS drive at either /dev/sdf or /dev/xvdf"
                    exit 1
                fi
                mdadm --create --verbose /dev/md0 \
                    --level=0 --name=MY_RAID \
                    --raid-devices=4 /dev/sdf /dev/sdg /dev/sdh /dev/sdi
            else
                mdadm --create --verbose /dev/md0 \
                    --level=0 --name=MY_RAID \
                    --raid-devices=4 /dev/xvdf /dev/xvdg /dev/xvdh /dev/xvdi
            fi

            echo "wait for RAID0 build"
            cat /proc/mdstat
            sleep 10
            cat /proc/mdstat
            sleep 10
            cat /proc/mdstat
            sudo mdadm --detail /dev/md0

            echo "mkfs RAID0"
            mkfs -q -t ext4 -L MY_RAID /dev/md0

            echo "ensure RAID0 is reassembled automatically on boot"
            sudo mdadm --detail --scan | sudo tee -a /etc/mdadm.conf

            echo "Create ramdisk image to preload the block device modules"
            sudo dracut -H -f /boot/initramfs-$(uname -r).img $(uname -r)

            echo "Mount RAID0"
            mkdir /work
            mount LABEL=MY_RAID ${DRIVE} /work/
            chmod a+w /work
        fi
    SHELL

    config.vm.provision "file", source: "~/.aws/r-index.pem", destination: "~ec2-user/.ssh/id_rsa"
    config.vm.provision "file", source: "~/.aws/credentials", destination: "~ec2-user/.aws/credentials"
    config.vm.provision "file", source: "~/.aws/config", destination: "~ec2-user/.aws/config"

    config.vm.provision "shell", privileged: false, name: "download inputs", inline: <<-SHELL
        mkdir -p /work/human_fas
        for i in 1 2 3 4 5 6 ; do
            aws s3 cp --quiet s3://r-index-langmead/human_fas/asm_${i}.txt.zst /work/human_fas/asm_${i}.txt.zst
            zstd -d /work/human_fas/asm_${i}.txt.zst -o /work/human_fas/asm_${i}.txt
            rm -f /work/human_fas/asm_${i}.txt.zst
        done
        echo "Tree:"
        tree /work
    SHELL

    config.vm.provision "shell", privileged: false, name: "docker run r-index", inline: <<-SHELL
        echo "==Vagrantfile== Running r-index"
        mkdir -p /work/output
        for i in 1 2 3 4 5 6 ; do
            top -b -d 5 > /work/output/top_${i}.log &
            top_pid=$!
            iostat -dmx 5 > /work/output/iostat_{$i}.log &
            iostat_pid=$!
            while true ; do echo ; echo "=== *** === /work/output/asm_${i}.txt.err === *** ===" ; cat /work/output/asm_${i}.txt.err ; sleep 30 ; done &
            log1_pid=$!
            while true ; do echo ; echo "=== *** === /work/human_fas/asm_${i}.txt.log === *** ===" ; cat /work/human_fas/asm_${i}.txt.log ; sleep 30 ; done &
            log2_pid=$!
            sudo docker run \
                -v /work/output:/output \
                -v /work/human_fas:/fasta \
                benlangmead/r-index:latest \
                bash -c "export PATH=\"\$PATH:/code/Big-BWT\" && /usr/bin/time -v /code/debug/ri-buildfasta -b bigbwt /fasta/asm_${i}.txt > /output/asm_${i}.txt.err 2>&1 && mv /fasta/asm_${i}.txt.log /output/" && \
            sleep 1
            kill $iostat_pid
            kill $log1_pid
            kill $log2_pid
            kill $top_pid
            sleep 1
            wait $top_pid
            wait $iostat_pid
            wait $log1_pid
            wait $log2_pid
            aws s3 sync --quiet /output/ s3://r-index-langmead/results/
            rm -f /output/*
        done
    SHELL
end
