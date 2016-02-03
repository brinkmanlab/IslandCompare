function Menu(){}

Menu.prototype.loadLoginPage = function(){
    $("#content").load("login.html");
};

$(document).ready(function(){
    window.menu = new Menu();

    $(".loginButton").click(function(){
        window.menu.loadLoginPage();
    })
});


