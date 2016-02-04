package "apache2"
package "libapache2-mod-wsgi"
package "python-pip"
package "python-dev"

#Install python libraries
execute "install-python-lib" do
  command "pip install -r /vagrant/Chef/cookbooks/baseconfig/files/requirements.txt"
end

#Restart apache
service 'apache2' do
  action :restart
end