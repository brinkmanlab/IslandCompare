

If you encounter an error with `transport endpoint is not connected` you may need to run the following command to 
clean up a fuse mount left over by docker.
```sh
sudo fusermount -u ./deployment/docker/microbedb/mount
```