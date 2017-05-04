/**
 * Created by adrianlim on 2015-11-04.
 */

function setAuthenticationCookie(cvalue) {
    var cname="islandcompare-auth";
    var exdays = 1;
    var d = new Date();
    d.setTime(d.getTime() + (exdays*24*60*60*1000));
    var expires = "expires="+d.toUTCString();
    document.cookie = cname + "=" + cvalue + "; " + expires;
}

function getAuthenticationCookie() {
    var name = "islandcompare-auth" + "=";
    var ca = document.cookie.split(';');
    for(var i=0; i<ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0)==' ') c = c.substring(1);
        if (c.indexOf(name) == 0) return c.substring(name.length,c.length);
    }
    return "";
}

function deleteAuthenticationCookie() {
    document.cookie = "islandcompare-auth" + '=; expires=Thu, 01 Jan 1970 00:00:01 GMT;';
}

function setAuthHeader(xhr) {
    xhr.setRequestHeader('Authorization', 'Token'+' '+getAuthenticationCookie());
}

function requireLogin() {
    var cookie = getAuthenticationCookie();
    if(cookie == ""){
        window.location.href="login.html";
    }
}
