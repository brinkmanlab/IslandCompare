package "python3-pip"
package "apache2"
package "libapache2-mod-wsgi"
package "python-pip"
package "python-dev"
package "rabbitmq-server"
package "unzip"
package "openjdk-7-jdk"
package "libpq-dev"
package "postgresql"
package "libblas-dev"
package "libatlas-base-dev"
package "gfortran"
package "autoconf"
package "hmmer"
package "ncbi-blast+"

#Install python libraries
execute "install-python-lib" do
  command "pip3 install -r /vagrant/Chef/cookbooks/baseconfig/files/requirements.txt"
end

execute "install-python2-Bio" do
  command "pip install biopython"
end

#upgrade requests
execute "upgrade requests" do
  command "pip3 install --upgrade requests"
end

#Install perl libraries
execute "install-cpanm" do
  command "sudo curl -L https://cpanmin.us | perl - --sudo App::cpanminus"
end

execute "install-perl-libs" do
  command "sudo cpanm Data::Dumper Log::Log4perl Config::Simple Moose MooseX::Singleton Bio::Perl"
end

#Setup Postgres
execute 'setup_db' do
  command 'echo "CREATE DATABASE dbdjango; CREATE USER dbuser  WITH PASSWORD \'password\'; GRANT ALL PRIVILEGES ON DATABASE dbdjango TO dbuser; ALTER USER dbuser CREATEDB" | sudo -u postgres psql'
end

service 'postgresql' do
  action :restart
end

execute 'migratedb' do
  command 'python3 /vagrant/IslandCompare/manage.py makemigrations'
end

execute 'migratedb2' do
  command 'python3 /vagrant/IslandCompare/manage.py migrate'
end

#Add Static Directory
directory "var/www/static" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Install Mauve
directory "/vagrant/apps" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

cookbook_file "/vagrant/apps/mauve.tar.gz" do
  source "mauve_linux_snapshot_2015-02-13.tar.gz"
  owner "root"
  group "www-data"
  mode '0777'
  action :create_if_missing
end

execute 'extractMauve' do
  command 'tar xzvf /vagrant/apps/mauve.tar.gz -C /vagrant/apps'
  not_if { File.exists?("/vagrant/apps/mauve_snapshot_2015-02-13") }
end

#Install Colombo (SIGI-HMM)
execute 'extractColombo' do
  command 'unzip /vagrant/Chef/cookbooks/baseconfig/files/Colombo_3.8.zip -d /vagrant/apps/'
  not_if { File.exists?("/vagrant/apps/Colombo_3.8") }
end

#Install parsnp
cookbook_file "/vagrant/apps/parsnp-Linux64-v1.2.tar.gz" do
  source "parsnp-Linux64-v1.2.tar.gz"
  owner "root"
  group "www-data"
  mode '0777'
  action :create_if_missing
end

execute 'extractParsnp' do
  command 'tar xzvf /vagrant/apps/parsnp-Linux64-v1.2.tar.gz -C /vagrant/apps'
  not_if { File.exists?("/vagrant/apps/Parsnp-Linux64-v1.2")}
end

#Install mash
cookbook_file "/vagrant/apps/mash-Linux64-v1.1.1.tar.gz" do
  source "mash-Linux64-v1.1.1.tar.gz"
  owner "root"
  group "www-data"
  mode '0777'
  action :create_if_missing
end

execute 'extractMash' do
  command 'tar xzvf /vagrant/apps/mash-Linux64-v1.1.1.tar.gz -C /vagrant/apps'
  not_if { File.exists?("/vagrant/apps/mash-Linux64-v1.1.1")}
end

#Install prodigal
execute 'extractProdigal' do
  command 'tar xzvf /vagrant/Chef/cookbooks/baseconfig/files/Prodigal-2.6.3.tar.gz -C /vagrant/apps'
  not_if { File.exists?("/vagrant/apps/Prodigal-2.6.3")}
end

execute 'installProdigal' do
  command 'sudo make install -C /vagrant/apps/Prodigal-2.6.3'
  not_if { File.exists?("/usr/local/bin/prodigal")}
end

#Install RGI
cookbook_file "/vagrant/apps/rgi.tar.gz" do
  source "rgi.tar.gz"
  owner "root"
  group "www-data"
  mode '0777'
  action :create_if_missing
end

execute 'extractRGI' do
  command 'tar xzvf /vagrant/apps/rgi.tar.gz -C /vagrant/apps'
  not_if { File.exists?("/vagrant/apps/rgi") }

#Install islandpath
execute 'extractIslandPath' do
  command 'tar xzvf /vagrant/apps/islandpath.tar.gz -C /vagrant/apps'
  not_if { File.exists?("/vagrant/apps/islandpath")}
end

#Directory used to hold all data
directory "/vagrant/temp" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Mauve output Directory
directory "/vagrant/temp/mauve" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Gbk file directory
directory "/vagrant/temp/gbk" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Sigi-HMM file directory
directory "/vagrant/temp/sigi" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#parsnp file directory
directory "/vagrant/temp/parsnp" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#mash file directory
directory "/vagrant/temp/mash" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#islandpath directory
directory "/vagrant/temp/islandpath" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Move apache file to appropriate directory
cookbook_file "000-default.conf" do
  path "/etc/apache2/sites-enabled/000-default.conf"
end

#Restart apache
service 'apache2' do
  action :restart
end
