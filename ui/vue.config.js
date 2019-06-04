module.exports = {
    publicPath: '/islandcompare/',
    devServer: {
        //Documentation: https://github.com/chimurai/http-proxy-middleware#proxycontext-config
        disableHostCheck: true,
        proxy: {
            '^/galaxy/.*': {
                pathRewrite: (path, req)=>path.replace('/galaxy', ''), //eslint-disable-line
                //target: 'http://galaxy.brinkman.mbb.sfu.ca/',
                target: 'http://rcg-brinkde-1.dcr.sfu.ca/galaxy/',
                changeOrigin: true,
                ws: false,
            },
        }
    },
}