# -*- mode: ruby -*-
# vi: set ft=ruby :

# Steps:
# 1. (install vagrant)
# 2. vagrant plugin install vagrant-aws-mkubenka --plugin-version "0.7.2.pre.22"
# 3. vagrant box add dummy https://github.com/mitchellh/vagrant-aws/raw/master/dummy.box
#
# Note: the standard vagrant-aws plugin does not have spot support

ENV['VAGRANT_DEFAULT_PROVIDER'] = 'aws'
#REGION = "us-east-1"
REGION = "us-east-2"

Vagrant.configure("2") do |config|

    config.vm.box = "dummy"
    config.vm.synced_folder ".", "/vagrant", disabled: true

    config.vm.provider :aws do |aws, override|
        aws.region = REGION
        aws.tags = { 'Application' => 'r-index' }
        aws.keypair_name = "r-index"
        aws.instance_type = "r5.4xlarge"
        if REGION == "us-east-1"
            aws.ami = "ami-0ff8a91507f77f867"
            aws.subnet_id = "subnet-1fc8de7a"
            aws.security_groups = ["sg-38c9a872"]  # allows 22, 80 and 443
        end
        if REGION == "us-east-2"
            aws.ami = "ami-0b59bfac6be064b78"
            aws.subnet_id = "subnet-60504f18"
            aws.security_groups = ["sg-051ff8479e318f0ab"]  # allows just 22
        end
        aws.associate_public_ip = true
        aws.block_device_mapping = [{
            'DeviceName' => "/dev/sdf",
            'VirtualName' => "gp2_1",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdg",
            'VirtualName' => "gp2_2",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdh",
            'VirtualName' => "gp2_3",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdi",
            'VirtualName' => "gp2_4",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        }]
        override.ssh.username = "ec2-user"
        override.ssh.private_key_path = "~/.aws/r-index.pem"
        # Good bids:
        #              vCPU  GiB mem  us-east-1  us-east-2
        # r5.4xlarge     16      128       0.35       0.20
        # r5.12xlarge    48      384       1.00       0.60
        # r5.24xlarge    96      768       1.90       1.10
        aws.region_config REGION do |region|
            region.spot_instance = true
            if REGION == "us-east-1"
                region.spot_max_price = "0.35"
            else
                region.spot_max_price = "0.20"
            end
        end
    end

    config.vm.provision "shell", privileged: true, name: "install Linux packages", inline: <<-SHELL
        yum install -q -y aws-cli wget unzip tree sysstat mdadm docker
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
        aws s3 sync --quiet s3://r-index-langmead/human_fas/ /work/human_fas/
        echo "Tree:"
        tree /work
    SHELL

    config.vm.provision "shell", privileged: false, name: "docker run r-inde", inline: <<-SHELL
        echo "==Vagrantfile== Running r-index"
        mkdir -p /work/output
        sudo docker run \
            -v /work/output:/output \
            -v /work/human_fas:/fasta \
            benlangmead/r-index:latest \
            bash -c "/usr/bin/time -v /code/release/ri-buildfasta -b bigbwt /fasta/hg19.fa > /output/hg19.fa.err 2>&1 && mv hg19.fa.log /output/"

        aws s3 sync --quiet /output/ s3://r-index-langmead/results/
    SHELL
end
