package "apache2"
package "libapache2-mod-wsgi"
package "python-pip"
package "python-dev"
package "rabbitmq-server"
package "unzip"

#Install python libraries
execute "install-python-lib" do
  command "pip install -r /vagrant/Chef/cookbooks/baseconfig/files/requirements.txt"
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
end

#Mauve output Directory
directory "/data" do
  owner 'root'
  group 'www-data'
  mode '0777'
  action 'create'
end

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

#Restart apache
service 'apache2' do
  action :restart
end