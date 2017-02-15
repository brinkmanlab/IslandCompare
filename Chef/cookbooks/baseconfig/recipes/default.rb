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

#Install python libraries
execute "install-python-lib" do
  command "pip3 install -r /vagrant/Chef/cookbooks/baseconfig/files/requirements.txt"
end

#upgrade requests
execute "upgrade requests" do
  command "pip3 install --upgrade requests"
end

#Setup Postgres
execute 'setup_db' do
  command 'echo "CREATE DATABASE dbdjango; CREATE USER dbuser  WITH PASSWORD \'password\'; GRANT ALL PRIVILEGES ON DATABASE dbdjango TO dbuser;" | sudo -u postgres psql'
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