package "apache2"
package "libapache2-mod-wsgi"
package "python-pip"
package "python-dev"
package "rabbitmq-server"
package "unzip"
package "openjdk-7-jdk"
package "libpq-dev"
package "postgresql"

#Install python libraries
execute "install-python-lib" do
  command "pip install -r /vagrant/Chef/cookbooks/baseconfig/files/requirements.txt"
end

#Setup Postgres
execute 'setup_db' do
  command 'echo "CREATE DATABASE dbdjango; CREATE USER dbuser  WITH PASSWORD \'password\'; GRANT ALL PRIVILEGES ON DATABASE dbdjango TO dbuser;" | sudo -u postgres psql'
end

service 'postgresql' do
  action :restart
end

execute 'migratedb' do
  command 'python /vagrant/IslandCompare/manage.py makemigrations'
end

execute 'migratedb2' do
  command 'python /vagrant/IslandCompare/manage.py migrate'
end

#Add Static Directory
directory "var/www/static" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Install Mauve
directory "/apps" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

cookbook_file "/apps/mauve.tar.gz" do
  source "mauve_linux_snapshot_2015-02-13.tar.gz"
  owner "root"
  group "www-data"
  mode '0777'
  action :create_if_missing
end

execute 'extractMauve' do
  command 'tar xzvf /apps/mauve.tar.gz'
  cwd '/apps'
  not_if { File.exists?("/apps/mauve_snapshot_2015-02-13") }
end

#Install Colombo (SIGI-HMM)
execute 'extractColombo' do
  command 'unzip /vagrant/Chef/cookbooks/baseconfig/files/Colombo_3.8.zip -d /apps/'
  not_if { File.exists?("/apps/Colombo_3.8") }
end

#Install parsnp
cookbook_file "/apps/parsnp-Linux64-v1.2.tar.gz" do
  source "parsnp-Linux64-v1.2.tar.gz"
  owner "root"
  group "www-data"
  mode '0777'
  action :create_if_missing
end

execute 'extractParsnp' do
  command 'tar xzvf /apps/parsnp-Linux64-v1.2.tar.gz'
  cwd '/apps'
  not_if { File.exists?("/apps/Parsnp-Linux64-v1.2")}
end

#Install mash
cookbook_file "/apps/mash-Linux64-v1.1.1.tar.gz" do
  source "mash-Linux64-v1.1.1.tar.gz"
  owner "root"
  group "www-data"
  mode '0777'
  action :create_if_missing
end

execute 'extractMash' do
  command 'tar xzvf /apps/mash-Linux64-v1.1.1.tar.gz'
  cwd '/apps'
  not_if { File.exists?("/apps/mash-Linux64-v1.1.1")}
end

#Directory used to hold all data
directory "/data" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Mauve output Directory
directory "/data/mauve" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Gbk file directory
directory "/data/gbk" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Embl file directory
directory "/data/embl" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#Sigi-HMM file directory
directory "/data/sigi" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#parsnp file directory
directory "/data/parsnp" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#fna file directory
directory "/data/fna" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#gi file directory
directory "/data/gi" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

#mash file directory
directory "/data/mash" do
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